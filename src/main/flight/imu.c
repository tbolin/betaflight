/*
 * This file is part of Cleanflight and Betaflight.
 *
 * Cleanflight and Betaflight are free software. You can redistribute
 * this software and/or modify this software under the terms of the
 * GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Cleanflight and Betaflight are distributed in the hope that they
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.
 *
 * If not, see <http://www.gnu.org/licenses/>.
 */

// Inertial Measurement Unit (IMU)

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "platform.h"

#include "build/build_config.h"
#include "build/debug.h"

#include "common/axis.h"
#include "common/vector.h"

#include "pg/pg.h"
#include "pg/pg_ids.h"

#include "drivers/time.h"

#include "fc/runtime_config.h"

#include "flight/gps_rescue.h"
#include "flight/imu.h"
#include "flight/mixer.h"
#include "flight/pid.h"

#include "io/gps.h"

#include "scheduler/scheduler.h"

#include "sensors/acceleration.h"
#include "sensors/barometer.h"
#include "sensors/compass.h"
#include "sensors/gyro.h"
#include "sensors/sensors.h"

#if defined(SIMULATOR_BUILD) && defined(SIMULATOR_MULTITHREAD)
#include <stdio.h>
#include <pthread.h>

static pthread_mutex_t imuUpdateLock;

#if defined(SIMULATOR_IMU_SYNC)
static uint32_t imuDeltaT = 0;
static bool imuUpdated = false;
#endif

#define IMU_LOCK pthread_mutex_lock(&imuUpdateLock)
#define IMU_UNLOCK pthread_mutex_unlock(&imuUpdateLock)

#else

#define IMU_LOCK
#define IMU_UNLOCK

#endif

// the limit (in degrees/second) beyond which we stop integrating
// omega_I. At larger spin rates the DCM PI controller can get 'dizzy'
// which results in false gyro drift. See
// https://drive.google.com/file/d/0ByvTkVQo3tqXQUVCVUNyZEgtRGs/view?usp=sharing&resourcekey=0-Mo4254cxdWWx2Y4mGN78Zw

#define SPIN_RATE_LIMIT 20

#define ATTITUDE_RESET_QUIET_TIME 250000   // 250ms - gyro quiet period after disarm before attitude reset
#define ATTITUDE_RESET_GYRO_LIMIT 15       // 15 deg/sec - gyro limit for quiet period
#define ATTITUDE_RESET_KP_GAIN    25.0     // dcmKpGain value to use during attitude reset
#define ATTITUDE_RESET_ACTIVE_TIME 500000  // 500ms - Time to wait for attitude to converge at high gain
#define GPS_COG_MIN_GROUNDSPEED 100        // 1.0m/s - the minimum groundspeed for a gps based IMU heading to be considered valid
                                           // Better to have some update than none for GPS Rescue at slow return speeds
#define GPS_COG_MAX_GROUNDSPEED 500        // 5.0m/s - Value for 'normal' 1.0 yaw IMU CogGain

bool canUseGPSHeading = true;

static float throttleAngleScale;
static int throttleAngleValue;
static float smallAngleCosZ = 0;

static imuRuntimeConfig_t imuRuntimeConfig;

float rMat[3][3];
static fpVector2_t north_ef;

#if defined(USE_ACC)
STATIC_UNIT_TESTED bool attitudeIsEstablished = false;
static float attitudeDeviation = 180.0f;
#endif

// quaternion of sensor frame relative to earth frame
STATIC_UNIT_TESTED quaternion q = QUATERNION_INITIALIZE;
STATIC_UNIT_TESTED quaternionProducts qP = QUATERNION_PRODUCTS_INITIALIZE;
// headfree quaternions
quaternion headfree = QUATERNION_INITIALIZE;
quaternion offset = QUATERNION_INITIALIZE;

// absolute angle inclination in multiple of 0.1 degree    180 deg = 1800
attitudeEulerAngles_t attitude = EULER_INITIALIZE;

PG_REGISTER_WITH_RESET_TEMPLATE(imuConfig_t, imuConfig, PG_IMU_CONFIG, 3);

#ifdef USE_RACE_PRO
#define DEFAULT_SMALL_ANGLE 180
#else
#define DEFAULT_SMALL_ANGLE 25
#endif

PG_RESET_TEMPLATE(imuConfig_t, imuConfig,
    .imu_dcm_kp = 2500,      // 1.0 * 10000
    .imu_dcm_ki = 0,         // 0.003 * 10000
    .small_angle = DEFAULT_SMALL_ANGLE,
    .imu_process_denom = 2,
    .mag_declination = 0
);

static void imuQuaternionComputeProducts(quaternion *quat, quaternionProducts *quatProd)
{
    quatProd->ww = quat->w * quat->w;
    quatProd->wx = quat->w * quat->x;
    quatProd->wy = quat->w * quat->y;
    quatProd->wz = quat->w * quat->z;
    quatProd->xx = quat->x * quat->x;
    quatProd->xy = quat->x * quat->y;
    quatProd->xz = quat->x * quat->z;
    quatProd->yy = quat->y * quat->y;
    quatProd->yz = quat->y * quat->z;
    quatProd->zz = quat->z * quat->z;
}

STATIC_UNIT_TESTED void imuComputeRotationMatrix(void)
{
    imuQuaternionComputeProducts(&q, &qP);

    rMat[0][0] = 1.0f - 2.0f * qP.yy - 2.0f * qP.zz;
    rMat[0][1] = 2.0f * (qP.xy + -qP.wz);
    rMat[0][2] = 2.0f * (qP.xz - -qP.wy);

    rMat[1][0] = 2.0f * (qP.xy - -qP.wz);
    rMat[1][1] = 1.0f - 2.0f * qP.xx - 2.0f * qP.zz;
    rMat[1][2] = 2.0f * (qP.yz + -qP.wx);

    rMat[2][0] = 2.0f * (qP.xz + -qP.wy);
    rMat[2][1] = 2.0f * (qP.yz - -qP.wx);
    rMat[2][2] = 1.0f - 2.0f * qP.xx - 2.0f * qP.yy;

#if defined(SIMULATOR_BUILD) && !defined(USE_IMU_CALC) && !defined(SET_IMU_FROM_EULER)
    rMat[1][0] = -2.0f * (qP.xy - -qP.wz);
    rMat[2][0] = -2.0f * (qP.xz + -qP.wy);
#endif
}

static float calculateThrottleAngleScale(uint16_t throttle_correction_angle)
{
    return (1800.0f / M_PIf) * (900.0f / throttle_correction_angle);
}

void imuConfigure(uint16_t throttle_correction_angle, uint8_t throttle_correction_value)
{
    imuRuntimeConfig.imuDcmKp = imuConfig()->imu_dcm_kp / 10000.0f;
    imuRuntimeConfig.imuDcmKi = imuConfig()->imu_dcm_ki / 10000.0f;
    // magnetic declination has negative sign (positive clockwise when seen from top)
    const float imuMagneticDeclinationRad = DEGREES_TO_RADIANS(imuConfig()->mag_declination / 10.0f);
    north_ef.x = cos_approx(imuMagneticDeclinationRad);
    north_ef.y = -sin_approx(imuMagneticDeclinationRad);

    smallAngleCosZ = cos_approx(degreesToRadians(imuConfig()->small_angle));

    throttleAngleScale = calculateThrottleAngleScale(throttle_correction_angle);

    throttleAngleValue = throttle_correction_value;
}

void imuInit(void)
{
#ifdef USE_GPS
    canUseGPSHeading = true;
#else
    canUseGPSHeading = false;
#endif

    imuComputeRotationMatrix();

#if defined(SIMULATOR_BUILD) && defined(SIMULATOR_MULTITHREAD)
    if (pthread_mutex_init(&imuUpdateLock, NULL) != 0) {
        printf("Create imuUpdateLock error!\n");
    }
#endif
}

#if defined(USE_ACC)
static float invSqrt(float x)
{
    return 1.0f / sqrtf(x);
}

STATIC_UNIT_TESTED void imuMahonyAHRSupdate(float dt, float gx, float gy, float gz,
                                bool useAcc, float ax, float ay, float az,
                                bool useMag,
                                float cogYawGain, float courseOverGround, const float dcmKpGain)
{
    static float integralFBx = 0.0f,  integralFBy = 0.0f, integralFBz = 0.0f;    // integral error terms scaled by Ki

    // Calculate general spin rate (rad/s)
    const float spin_rate = sqrtf(sq(gx) + sq(gy) + sq(gz));

    // Use raw heading error (from GPS or whatever else)
    float ex = 0, ey = 0, ez = 0;
    if (cogYawGain != 0.0f) {
        // Used in a GPS Rescue to boost IMU yaw gain when course over ground and velocity to home differ significantly

        // Compute heading vector in EF from scalar CoG. CoG is clockwise from North
        // Note that Earth frame X is pointing north and sin/cos argument is anticlockwise
        const fpVector2_t cog_ef = {.x = cos_approx(-courseOverGround), .y = sin_approx(-courseOverGround)};
#define THRUST_COG 1
#if THRUST_COG
        const fpVector2_t heading_ef = {.x = rMat[X][Z], .y = rMat[Y][Z]};  // body Z axis (up) - direction of thrust vector
#else
        const fpVector2_t heading_ef = {.x = rMat[0][0], .y = rMat[1][0]};  // body X axis. Projected vector magnitude is reduced as pitch increases
#endif
        // cross product = 1 * |heading| * sin(angle) (magnitude of Z vector in 3D)
        // order operands so that rotation is in direction of zero error
        const float cross = vector2Cross(&heading_ef, &cog_ef);
        // dot product, 1 * |heading| * cos(angle)
        const float dot = vector2Dot(&heading_ef, &cog_ef);
        // use cross product / sin(angle) when error < 90deg (cos > 0),
        //   |heading| if error is larger (cos < 0)
        const float heading_mag = vector2Mag(&heading_ef);
        float ez_ef = (dot > 0) ? cross : (cross < 0 ? -1.0f : 1.0f) * heading_mag;
#if THRUST_COG
        // increase gain for small tilt (just heuristic; sqrt is cheap on F4+)
        ez_ef /= sqrtf(heading_mag);
#endif
        ez_ef *= cogYawGain;          // apply gain parameter
        // convert to body frame
        ex += rMat[2][0] * ez_ef;
        ey += rMat[2][1] * ez_ef;
        ez += rMat[2][2] * ez_ef;

        DEBUG_SET(DEBUG_ATTITUDE, 3, (ez_ef * 100));
    }

    DEBUG_SET(DEBUG_ATTITUDE, 2, cogYawGain * 100.0f);
    DEBUG_SET(DEBUG_ATTITUDE, 7, (dcmKpGain * 100));

#ifdef USE_MAG
    // Use measured magnetic field vector
    fpVector3_t mag_bf = {{mag.magADC[X], mag.magADC[Y], mag.magADC[Z]}};
    float magNormSquared = vectorNormSquared(&mag_bf);
    fpVector3_t mag_ef;
    matrixVectorMul(&mag_ef, (const fpMat33_t*)&rMat, &mag_bf); // BF->EF true north

#ifdef USE_GPS_RESCUE
    // Encapsulate additional operations in a block so that it is only executed when the according debug mode is used
    // Only re-calculate magYaw when there is a new Mag data reading, to avoid spikes
    if (debugMode == DEBUG_GPS_RESCUE_HEADING && mag.isNewMagADCFlag) {
        fpMat33_t rMatZTrans;
        yawToRotationMatrixZ(&rMatZTrans, -atan2_approx(rMat[1][0], rMat[0][0]));
        fpVector3_t mag_ef_yawed;
        matrixVectorMul(&mag_ef_yawed, &rMatZTrans, &mag_ef); // EF->EF yawed
        // Magnetic yaw is the angle between true north and the X axis of the body frame
        int16_t magYaw = lrintf((atan2_approx(mag_ef_yawed.y, mag_ef_yawed.x) * (1800.0f / M_PIf)));
        if (magYaw < 0) {
            magYaw += 3600;
        }
        DEBUG_SET(DEBUG_GPS_RESCUE_HEADING, 4, magYaw); // mag heading in degrees * 10
        // reset new mag data flag to false to initiate monitoring for new Mag data.
        // note that if the debug doesn't run, this reset will not occur, and we won't waste cycles on the comparison
        mag.isNewMagADCFlag = false;
    }
#endif

    if (useMag && magNormSquared > 0.01f) {
        // Normalise magnetometer measurement
        vectorNormalize(&mag_ef, &mag_ef);

        // For magnetometer correction we make an assumption that magnetic field is perpendicular to gravity (ignore Z-component in EF).
        // This way magnetic field will only affect heading and wont mess roll/pitch angles

        // (hx; hy; 0) - measured mag field vector in EF (forcing Z-component to zero)
        // (bx; 0; 0) - reference mag field vector heading due North in EF (assuming Z-component is zero)
        mag_ef.z = 0.0f;                // project to XY plane (optimized away)

        // magnetometer error is cross product between estimated magnetic north and measured magnetic north (calculated in EF)
        // increase gain on large misalignment
        const float dot = vector2Dot((fpVector2_t*)&mag_ef, &north_ef);
        const float cross = vector2Cross((fpVector2_t*)&mag_ef, &north_ef);
        const float ez_ef = (dot > 0) ? cross : (cross < 0 ? -1.0f : 1.0f) * vector2Mag((fpVector2_t*)&mag_ef);
        // Rotate mag error vector back to BF and accumulate
        ex += rMat[2][0] * ez_ef * dcmKpGain;
        ey += rMat[2][1] * ez_ef * dcmKpGain;
        ez += rMat[2][2] * ez_ef * dcmKpGain;
    }
#else
    UNUSED(useMag);
#endif

    // Use measured acceleration vector
    float recipAccNorm = sq(ax) + sq(ay) + sq(az);
    if (useAcc && recipAccNorm > 0.01f) {
        // Normalise accelerometer measurement; useAcc is true when all smoothed acc axes are within 20% of 1G
        recipAccNorm = invSqrt(recipAccNorm);
        ax *= recipAccNorm;
        ay *= recipAccNorm;
        az *= recipAccNorm;

        fpVector3_t acc_upp = {.x = ax, .y = ay, .z = az};
        fpVector3_t est_upp = {.x = rMat[2][0], .y = rMat[2][1], .z = rMat[2][2]};

        // Error is sum of cross product between estimated direction and measured direction of gravity
        fpVector3_t acc_err;
        vectorCrossProduct(&acc_err, &acc_upp, &est_upp);

        const float dot = vector3Dot(&acc_upp, &est_upp);

        // To avoid the gain decreasing for angles > 90 degrees:
        // set magnitude of error vector to 1 + |cos| of angle between estimated
        // and measured downwards vectors if the absolute angle of the error is > 90 degrees
        acc_err = (dot > 0) ? acc_err : *vectorScale(&acc_err, vectorNormalize(&acc_err, &acc_err), 1.0f - dot);

        ex += acc_err.x * dcmKpGain;
        ey += acc_err.y * dcmKpGain;
        ez += acc_err.z * dcmKpGain;
    }

    // Compute and apply integral feedback if enabled
    if (imuRuntimeConfig.imuDcmKi > 0.0f) {
        // Stop integrating if spinning beyond the certain limit
        if (spin_rate < DEGREES_TO_RADIANS(SPIN_RATE_LIMIT)) {
            const float dcmKiGain = imuRuntimeConfig.imuDcmKi;
            integralFBx += dcmKiGain * ex * dt;    // integral error scaled by Ki
            integralFBy += dcmKiGain * ey * dt;
            integralFBz += dcmKiGain * ez * dt;
        }
    } else {
        integralFBx = 0.0f;    // prevent integral windup
        integralFBy = 0.0f;
        integralFBz = 0.0f;
    }

    // Apply proportional and integral feedback
    gx += ex + integralFBx;
    gy += ey + integralFBy;
    gz += ez + integralFBz;

#define IMU_USE_4TH_ORDER_INTEGRATION 0
#if !IMU_USE_4TH_ORDER_INTEGRATION
    // Integrate rate of change of quaternion
    gx *= (0.5f * dt);
    gy *= (0.5f * dt);
    gz *= (0.5f * dt);

    quaternion buffer;
    buffer.w = q.w;
    buffer.x = q.x;
    buffer.y = q.y;
    buffer.z = q.z;

    q.w += (-buffer.x * gx - buffer.y * gy - buffer.z * gz);
    q.x += (+buffer.w * gx + buffer.y * gz - buffer.z * gy);
    q.y += (+buffer.w * gy - buffer.x * gz + buffer.z * gx);
    q.z += (+buffer.w * gz + buffer.x * gy - buffer.y * gx);

#else
    const quaternion k1 = {
        .w = (-q.x * gx - q.y * gy - q.z * gz) * 0.5f,
        .x = (+q.w * gx + q.y * gz - q.z * gy) * 0.5f,
        .y = (+q.w * gy - q.x * gz + q.z * gx) * 0.5f,
        .z = (+q.w * gz + q.x * gy - q.y * gx) * 0.5f,
    };

    const quaternion b1 = {
        .w = q.w + k1.w * 0.5f * dt,
        .x = q.x + k1.x * 0.5f * dt,
        .y = q.y + k1.y * 0.5f * dt,
        .z = q.z + k1.z * 0.5f * dt,
    };

    const quaternion k2 = {
        .w = (-b1.x * gx - b1.y * gy - b1.z * gz) * 0.5f,
        .x = (+b1.w * gx + b1.y * gz - b1.z * gy) * 0.5f,
        .y = (+b1.w * gy - b1.x * gz + b1.z * gx) * 0.5f,
        .z = (+b1.w * gz + b1.x * gy - b1.y * gx) * 0.5f,
    };

    const quaternion b2 = {
        .w = q.w + k2.w * 0.5 * dt,
        .x = q.x + k2.x * 0.5 * dt,
        .y = q.y + k2.y * 0.5 * dt,
        .z = q.z + k2.z * 0.5 * dt,
    };

    const quaternion k3 = {
        .w = (-b2.x * gx - b2.y * gy - b2.z * gz) * 0.5f,
        .x = (+b2.w * gx + b2.y * gz - b2.z * gy) * 0.5f,
        .y = (+b2.w * gy - b2.x * gz + b2.z * gx) * 0.5f,
        .z = (+b2.w * gz + b2.x * gy - b2.y * gx) * 0.5f,
    };

    const quaternion b3 = {
        .w = q.w + k3.w * dt,
        .x = q.x + k3.x * dt,
        .y = q.y + k3.y * dt,
        .z = q.z + k3.z * dt,
    };

    const quaternion k4 = {
        .w = (-b3.x * gx - b3.y * gy - b3.z * gz) * 0.5f,
        .x = (+b3.w * gx + b3.y * gz - b3.z * gy) * 0.5f,
        .y = (+b3.w * gy - b3.x * gz + b3.z * gx) * 0.5f,
        .z = (+b3.w * gz + b3.x * gy - b3.y * gx) * 0.5f,
    };

    q.w += (k1.w + 2.0f * k2.w + 2.0f * k3.w + k4.w) * dt / 6.0f;
    q.x += (k1.x + 2.0f * k2.x + 2.0f * k3.x + k4.x) * dt / 6.0f;
    q.y += (k1.y + 2.0f * k2.y + 2.0f * k3.y + k4.y) * dt / 6.0f;
    q.z += (k1.z + 2.0f * k2.z + 2.0f * k3.z + k4.z) * dt / 6.0f;
    DEBUG_SET(DEBUG_IMU_GAIN, 2, lrintf(10000.0f * sqrtf(sq(q.w - b1.w) + sq(q.x - b1.x) + sq(q.y - b1.y) + sq(q.z - b1.z))));
#endif

    // Normalise quaternion
    float recipNorm = invSqrt(sq(q.w) + sq(q.x) + sq(q.y) + sq(q.z));
    q.w *= recipNorm;
    q.x *= recipNorm;
    q.y *= recipNorm;
    q.z *= recipNorm;

    // Pre-compute rotation matrix from quaternion
    imuComputeRotationMatrix();

    attitudeIsEstablished = true;
}

STATIC_UNIT_TESTED void imuUpdateEulerAngles(void)
{
    quaternionProducts buffer;

    if (FLIGHT_MODE(HEADFREE_MODE)) {
       imuQuaternionComputeProducts(&headfree, &buffer);

       attitude.values.roll = lrintf(atan2_approx((+2.0f * (buffer.wx + buffer.yz)), (+1.0f - 2.0f * (buffer.xx + buffer.yy))) * (1800.0f / M_PIf));
       attitude.values.pitch = lrintf(((0.5f * M_PIf) - acos_approx(+2.0f * (buffer.wy - buffer.xz))) * (1800.0f / M_PIf));
       attitude.values.yaw = lrintf((-atan2_approx((+2.0f * (buffer.wz + buffer.xy)), (+1.0f - 2.0f * (buffer.yy + buffer.zz))) * (1800.0f / M_PIf)));
    } else {
       attitude.values.roll = lrintf(atan2_approx(rMat[2][1], rMat[2][2]) * (1800.0f / M_PIf));
       attitude.values.pitch = lrintf(((0.5f * M_PIf) - acos_approx(-rMat[2][0])) * (1800.0f / M_PIf));
       attitude.values.yaw = lrintf((-atan2_approx(rMat[1][0], rMat[0][0]) * (1800.0f / M_PIf)));
    }

    if (attitude.values.yaw < 0) {
        attitude.values.yaw += 3600;
    }
}

static float imuUpdateDeviation(const float deviation, const float imuDt, float durationSaturated)
{
    // Increase estiumated deviation based on time spent integrating
    static const float err_rate_saturated = 500.0f; // guestimated gyro drift in deg/s when saturated
    static const float err_rate_normal = 1.0f;  // static gyro drift in deg/s under normal circumstances
    if (durationSaturated > imuDt) { durationSaturated = imuDt; }
    const float normal_duration = imuDt - durationSaturated;
    const float accumulated_error = err_rate_normal * normal_duration + err_rate_saturated * durationSaturated;
    return constrainf(deviation + accumulated_error, 0.0f, 180.0f);
}

static float imuAccTrust(const float* accAverage, const float* gyroAverage, const float gRecip)
{
    const fpVector3_t accVector = {.x = accAverage[X], .y = accAverage[Y], .z = accAverage[Z]};
    const fpVector3_t gyroVector = {.x = gyroAverage[X], .y = gyroAverage[Y], .z = gyroAverage[Z]};
    const float accNorm = vectorNorm(&accVector) * gRecip;
    const float gyroNorm = vectorNorm(&gyroVector);
    DEBUG_SET(DEBUG_IMU_GAIN, 4, lrintf(accNorm * 100.0f));
    DEBUG_SET(DEBUG_IMU_GAIN, 5, lrintf(gyroNorm * 1.0f));
    return tent(accNorm, 1.0f, 0.2f) * tent(gyroNorm, 0.0f, 100.0f);
}

// Calculate gain based on
// * how close the magnitude of the acc is to 1g
// * how close the rotational rate is to 0
// * how much we trust the current estimate
static float imuCalcKpGain(float *deviation, const float baseKp, const float dt, const float accTrust)
{
    const float deviationGainMax = 100.0f * baseKp;
    const float deviationGainMin = 0.1f * baseKp;
    const float kpGain = scaleRangef(*deviation, 0.0f, 180.0f, deviationGainMin, deviationGainMax) * accTrust;

    // decrease the estimated deviation based on how much we currently trust the acc
    static const float heuristic = 1.0f;
    const float omega = 2.0f * M_PIf * dt * (heuristic * baseKp * accTrust);
    *deviation *= 1.0f - omega / (omega + 1.0f);
    DEBUG_SET(DEBUG_IMU_GAIN, 0, lrintf(kpGain * 100.0f));
    DEBUG_SET(DEBUG_IMU_GAIN, 1, lrintf(*deviation));
    // DEBUG_SET(DEBUG_IMU_GAIN, 2, lrintf(dt * 100000.0f));
    DEBUG_SET(DEBUG_IMU_GAIN, 3, lrintf(accTrust * 1000.0f));
    return kpGain;
}

#if defined(USE_GPS)
static void imuComputeQuaternionFromRPY(quaternionProducts *quatProd, int16_t initialRoll, int16_t initialPitch, int16_t initialYaw)
{
    if (initialRoll > 1800) {
        initialRoll -= 3600;
    }

    if (initialPitch > 1800) {
        initialPitch -= 3600;
    }

    if (initialYaw > 1800) {
        initialYaw -= 3600;
    }

    const float cosRoll = cos_approx(DECIDEGREES_TO_RADIANS(initialRoll) * 0.5f);
    const float sinRoll = sin_approx(DECIDEGREES_TO_RADIANS(initialRoll) * 0.5f);

    const float cosPitch = cos_approx(DECIDEGREES_TO_RADIANS(initialPitch) * 0.5f);
    const float sinPitch = sin_approx(DECIDEGREES_TO_RADIANS(initialPitch) * 0.5f);

    const float cosYaw = cos_approx(DECIDEGREES_TO_RADIANS(-initialYaw) * 0.5f);
    const float sinYaw = sin_approx(DECIDEGREES_TO_RADIANS(-initialYaw) * 0.5f);

    const float q0 = cosRoll * cosPitch * cosYaw + sinRoll * sinPitch * sinYaw;
    const float q1 = sinRoll * cosPitch * cosYaw - cosRoll * sinPitch * sinYaw;
    const float q2 = cosRoll * sinPitch * cosYaw + sinRoll * cosPitch * sinYaw;
    const float q3 = cosRoll * cosPitch * sinYaw - sinRoll * sinPitch * cosYaw;

    quatProd->xx = sq(q1);
    quatProd->yy = sq(q2);
    quatProd->zz = sq(q3);

    quatProd->xy = q1 * q2;
    quatProd->xz = q1 * q3;
    quatProd->yz = q2 * q3;

    quatProd->wx = q0 * q1;
    quatProd->wy = q0 * q2;
    quatProd->wz = q0 * q3;

    imuComputeRotationMatrix();

    attitudeIsEstablished = true;
}
#endif

static void imuCalculateEstimatedAttitude(timeUs_t currentTimeUs)
{
    static timeUs_t previousIMUUpdateTime;
    bool useMag = false;
    float cogYawGain = 0.0f; // IMU yaw gain to be applied in imuMahonyAHRSupdate from ground course, default to no correction from CoG
    float courseOverGround = 0; // To be used when cogYawGain is non-zero, in radians

    const timeDelta_t deltaT = currentTimeUs - previousIMUUpdateTime;
    previousIMUUpdateTime = currentTimeUs;

#ifdef USE_MAG
    if (sensors(SENSOR_MAG) && compassIsHealthy()
#ifdef USE_GPS_RESCUE
        && !gpsRescueDisableMag()
#endif
        ) {
        useMag = true;
    }
#endif
#if defined(USE_GPS)
    if (!useMag && sensors(SENSOR_GPS) && STATE(GPS_FIX) && gpsSol.numSat > GPS_MIN_SAT_COUNT && gpsSol.groundSpeed >= GPS_COG_MIN_GROUNDSPEED) {
        // Use GPS course over ground to correct attitude.values.yaw
        const float imuYawGroundspeed = fminf (gpsSol.groundSpeed / GPS_COG_MAX_GROUNDSPEED, 2.0f);
        courseOverGround = DECIDEGREES_TO_RADIANS(gpsSol.groundCourse);
        cogYawGain = (FLIGHT_MODE(GPS_RESCUE_MODE)) ? gpsRescueGetImuYawCogGain() : imuYawGroundspeed;
        cogYawGain *= imuRuntimeConfig.dcm_kp;
        // normally update yaw heading with GPS data, but when in a Rescue, modify the IMU yaw gain dynamically
        if (shouldInitializeGPSHeading()) {
            // Reset our reference and reinitialize quaternion.
            // shouldInitializeGPSHeading() returns true only once.
            imuComputeQuaternionFromRPY(&qP, attitude.values.roll, attitude.values.pitch, gpsSol.groundCourse);
            cogYawGain = 0.0f; // Don't use the COG when we first initialize
        }
    }
#endif

#if defined(SIMULATOR_BUILD) && !defined(USE_IMU_CALC)
    UNUSED(imuMahonyAHRSupdate);
    UNUSED(attitudeDeviation);
    UNUSED(useMag);
    UNUSED(cogYawGain);
    UNUSED(canUseGPSHeading);
    UNUSED(courseOverGround);
    UNUSED(deltaT);
    UNUSED(imuAccTrust);
    UNUSED(imuCalcKpGain);
    UNUSED(imuUpdateDeviation);
#else

#if defined(SIMULATOR_BUILD) && defined(SIMULATOR_IMU_SYNC)
//  printf("[imu]deltaT = %u, imuDeltaT = %u, currentTimeUs = %u, micros64_real = %lu\n", deltaT, imuDeltaT, currentTimeUs, micros64_real());
    deltaT = imuDeltaT;
#endif
    float gyroAverage[XYZ_AXIS_COUNT];
    for (int axis = 0; axis < XYZ_AXIS_COUNT; ++axis) {
        gyroAverage[axis] = gyroGetFilteredDownsampled(axis);
    }

    const float dt = deltaT * 1e-6f;

    attitudeDeviation = imuUpdateDeviation(attitudeDeviation, dt, gyroGetDurationSpentSaturated());

    const float accTrust = imuAccTrust(acc.accADC, gyroAverage, acc.dev.acc_1G_rec);
    const bool useAcc = accTrust > 0.0f;

    const float adjustedAccKp = imuCalcKpGain(&attitudeDeviation, imuRuntimeConfig.dcm_kp, dt, accTrust);

    imuMahonyAHRSupdate(dt,
                        DEGREES_TO_RADIANS(gyroAverage[X]), DEGREES_TO_RADIANS(gyroAverage[Y]), DEGREES_TO_RADIANS(gyroAverage[Z]),
                        useAcc, acc.accADC[X], acc.accADC[Y], acc.accADC[Z],
                        useMag,
                        cogYawGain, courseOverGround,  adjustedAccKp);

    imuUpdateEulerAngles();
#endif
}

static int calculateThrottleAngleCorrection(void)
{
    /*
    * Use 0 as the throttle angle correction if we are inverted, vertical or with a
    * small angle < 0.86 deg
    * TODO: Define this small angle in config.
    */
    if (getCosTiltAngle() <= 0.015f) {
        return 0;
    }
    int angle = lrintf(acos_approx(getCosTiltAngle()) * throttleAngleScale);
    if (angle > 900)
        angle = 900;
    return lrintf(throttleAngleValue * sin_approx(angle / (900.0f * M_PIf / 2.0f)));
}

void imuUpdateAttitude(timeUs_t currentTimeUs)
{
    if (sensors(SENSOR_ACC) && acc.isAccelUpdatedAtLeastOnce && gyroIsCalibrationComplete()) {
        IMU_LOCK;
#if defined(SIMULATOR_BUILD) && defined(SIMULATOR_IMU_SYNC)
        if (imuUpdated == false) {
            IMU_UNLOCK;
            return;
        }
        imuUpdated = false;
#endif
        imuCalculateEstimatedAttitude(currentTimeUs);
        IMU_UNLOCK;

        // Update the throttle correction for angle and supply it to the mixer
        int throttleAngleCorrection = 0;
        if (throttleAngleValue && (FLIGHT_MODE(ANGLE_MODE) || FLIGHT_MODE(HORIZON_MODE)) && ARMING_FLAG(ARMED)) {
            throttleAngleCorrection = calculateThrottleAngleCorrection();
        }
        mixerSetThrottleAngleCorrection(throttleAngleCorrection);

    } else {
        acc.accADC[X] = 0;
        acc.accADC[Y] = 0;
        acc.accADC[Z] = 0;
        schedulerIgnoreTaskStateTime();
    }

    DEBUG_SET(DEBUG_ATTITUDE, X, acc.accADC[X]); // roll
    DEBUG_SET(DEBUG_ATTITUDE, Y, acc.accADC[Y]); // pitch
}
#endif // USE_ACC

bool shouldInitializeGPSHeading(void)
{
    static bool initialized = false;

    if (!initialized) {
        initialized = true;

        return true;
    }

    return false;
}

float getCosTiltAngle(void)
{
    return rMat[2][2];
}

void getQuaternion(quaternion *quat)
{
   quat->w = q.w;
   quat->x = q.x;
   quat->y = q.y;
   quat->z = q.z;
}

#ifdef SIMULATOR_BUILD
void imuSetAttitudeRPY(float roll, float pitch, float yaw)
{
    IMU_LOCK;

    attitude.values.roll = roll * 10;
    attitude.values.pitch = pitch * 10;
    attitude.values.yaw = yaw * 10;

    IMU_UNLOCK;
}

void imuSetAttitudeQuat(float w, float x, float y, float z)
{
    IMU_LOCK;

    q.w = w;
    q.x = x;
    q.y = y;
    q.z = z;

    imuComputeRotationMatrix();

    attitudeIsEstablished = true;

    imuUpdateEulerAngles();

    IMU_UNLOCK;
}
#endif
#if defined(SIMULATOR_BUILD) && defined(SIMULATOR_IMU_SYNC)
void imuSetHasNewData(uint32_t dt)
{
    IMU_LOCK;

    imuUpdated = true;
    imuDeltaT = dt;

    IMU_UNLOCK;
}
#endif

bool imuQuaternionHeadfreeOffsetSet(void)
{
    if ((abs(attitude.values.roll) < 450)  && (abs(attitude.values.pitch) < 450)) {
        const float yaw = -atan2_approx((+2.0f * (qP.wz + qP.xy)), (+1.0f - 2.0f * (qP.yy + qP.zz)));

        offset.w = cos_approx(yaw/2);
        offset.x = 0;
        offset.y = 0;
        offset.z = sin_approx(yaw/2);

        return true;
    } else {
        return false;
    }
}

void imuQuaternionMultiplication(quaternion *q1, quaternion *q2, quaternion *result)
{
    const float A = (q1->w + q1->x) * (q2->w + q2->x);
    const float B = (q1->z - q1->y) * (q2->y - q2->z);
    const float C = (q1->w - q1->x) * (q2->y + q2->z);
    const float D = (q1->y + q1->z) * (q2->w - q2->x);
    const float E = (q1->x + q1->z) * (q2->x + q2->y);
    const float F = (q1->x - q1->z) * (q2->x - q2->y);
    const float G = (q1->w + q1->y) * (q2->w - q2->z);
    const float H = (q1->w - q1->y) * (q2->w + q2->z);

    result->w = B + (- E - F + G + H) / 2.0f;
    result->x = A - (+ E + F + G + H) / 2.0f;
    result->y = C + (+ E - F + G - H) / 2.0f;
    result->z = D + (+ E - F - G + H) / 2.0f;
}

void imuQuaternionHeadfreeTransformVectorEarthToBody(t_fp_vector_def *v)
{
    quaternionProducts buffer;

    imuQuaternionMultiplication(&offset, &q, &headfree);
    imuQuaternionComputeProducts(&headfree, &buffer);

    const float x = (buffer.ww + buffer.xx - buffer.yy - buffer.zz) * v->X + 2.0f * (buffer.xy + buffer.wz) * v->Y + 2.0f * (buffer.xz - buffer.wy) * v->Z;
    const float y = 2.0f * (buffer.xy - buffer.wz) * v->X + (buffer.ww - buffer.xx + buffer.yy - buffer.zz) * v->Y + 2.0f * (buffer.yz + buffer.wx) * v->Z;
    const float z = 2.0f * (buffer.xz + buffer.wy) * v->X + 2.0f * (buffer.yz - buffer.wx) * v->Y + (buffer.ww - buffer.xx - buffer.yy + buffer.zz) * v->Z;

    v->X = x;
    v->Y = y;
    v->Z = z;
}

bool isUpright(void)
{
#ifdef USE_ACC
    return !sensors(SENSOR_ACC) || (attitudeIsEstablished && getCosTiltAngle() > smallAngleCosZ);
#else
    return true;
#endif
}
