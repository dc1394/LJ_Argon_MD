/*! \file Ar_moleculardynamics.h
	\brief アルゴンに対して、分子動力学シミュレーションを行うクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#include "DXUT.h"
#include "Ar_moleculardynamics.h"
#include "../myrandom/myrand.h"
#include <tbb/parallel_for.h>           // for tbb::parallel_for
#include <tbb/partitioner.h>            // for tbb::auto_partitioner
#include <tbb/task_scheduler_init.h>    // for tbb::task_scheduler_init

namespace moleculardynamics {
    // #region static private 定数

	double const Ar_moleculardynamics::FIRSTSCALE = 1.0;

	double const Ar_moleculardynamics::FIRSTTEMP = 50.0;

	double const Ar_moleculardynamics::ALPHA = 0.2;

	double const Ar_moleculardynamics::AVOGADRO_CONSTANT = 6.022140857E+23;

	double const Ar_moleculardynamics::DT = 0.001;

    double const Ar_moleculardynamics::KB = 1.3806488E-23;

    double const Ar_moleculardynamics::SIGMA = 3.405E-10;

	double const Ar_moleculardynamics::TAU =
		std::sqrt(0.039948 / Ar_moleculardynamics::AVOGADRO_CONSTANT * Ar_moleculardynamics::SIGMA * Ar_moleculardynamics::SIGMA / Ar_moleculardynamics::YPSILON);

    double const Ar_moleculardynamics::YPSILON = 1.6540172624E-21;

    // #endregion static private 定数

	// #region コンストラクタ

	Ar_moleculardynamics::Ar_moleculardynamics()
		:
		lat([this] { return lat_; }, nullptr),
		MD_iter([this] { return MD_iter_; }, nullptr),
		X([this] { return std::cref(X_); }, nullptr),
		Y([this] { return std::cref(Y_); }, nullptr),
		Z([this] { return std::cref(Z_); }, nullptr),
		dt2(DT * DT),
		FX(Nc * Nc * Nc * 4),
		FY(Nc * Nc * Nc * 4),
		FZ(Nc * Nc * Nc * 4),
		rc2(rc * rc),
		rcm6(std::pow(rc, -6.0)),
		rcm12(std::pow(rc, -12.0)),
		Tg_(Ar_moleculardynamics::FIRSTTEMP * Ar_moleculardynamics::KB / Ar_moleculardynamics::YPSILON),
		Vrc(4.0 * (rcm12 - rcm6)),
		VX(Nc * Nc * Nc * 4),
		VY(Nc * Nc * Nc * 4),
		VZ(Nc * Nc * Nc * 4),
		X_(Nc * Nc * Nc * 4),
		X1(Nc * Nc * Nc * 4),
		Y_(Nc * Nc * Nc * 4),
		Y1(Nc * Nc * Nc * 4),
		Z_(Nc * Nc * Nc * 4),
		Z1(Nc * Nc * Nc * 4)
	{
		// initalize parameters
		lat_ = std::pow(2.0, 2.0 / 3.0) * scale_;

        recalc();

		periodiclen = lat_ * static_cast<double>(Nc);
		
		tbb::task_scheduler_init init;
	}

	// #endregion コンストラクタ

    // #region publicメンバ関数

    void Ar_moleculardynamics::Calc_Forces()
    {
        for (int n = 0; n < NumAtom; n++) {
            FX[n] = 0.0;
            FY[n] = 0.0;
            FZ[n] = 0.0;
        }

		tbb::parallel_for(
			0,
			NumAtom,
			1,
			[this](std::int32_t n) {
			for (auto m = 0; m < NumAtom; m++) {

				// ±ncp分のセル内の原子との相互作用を計算
				for (auto i = -ncp; i <= ncp; i++) {
					for (auto j = -ncp; j <= ncp; j++) {
						for (auto k = -ncp; k <= ncp; k++) {
							auto const sx = static_cast<double>(i) * periodiclen;
							auto const sy = static_cast<double>(j) * periodiclen;
							auto const sz = static_cast<double>(k) * periodiclen;

							// 自分自身との相互作用を排除
							if (n != m || i != 0 || j != 0 || k != 0) {
								auto const dx = X_[n] - (X_[m] + sx);
								auto const dy = Y_[n] - (Y_[m] + sy);
								auto const dz = Z_[n] - (Z_[m] + sz);

								auto const r2 = norm2(dx, dy, dz);
								// 打ち切り距離内であれば計算
								if (r2 <= rc2) {
									auto const r = sqrt(r2);
									auto const rm6 = 1.0 / (r2 * r2 * r2);
									auto const rm7 = rm6 / r;
									auto const rm12 = rm6 * rm6;
									auto const rm13 = rm12 / r;

									auto const Fr = 48.0 * rm13 - 24.0 * rm7;

									FX[n] += dx / r * Fr;
									FY[n] += dy / r * Fr;
									FZ[n] += dz / r * Fr;
								}
							}
						}
					}
				}
			}
		},
		tbb::auto_partitioner());
    }
	
	float Ar_moleculardynamics::getDeltat() const
	{
		return static_cast<float>(t * TAU * 1.0E+12);
	}

    float Ar_moleculardynamics::getForce(std::int32_t n) const
    {
        return static_cast<float>(std::sqrt(norm2(FX[n], FY[n], FZ[n])));
    }

    double Ar_moleculardynamics::getTcalc() const
    {
        return Tc_ * Ar_moleculardynamics::YPSILON / Ar_moleculardynamics::KB;
    }

    double Ar_moleculardynamics::getTgiven() const
    {
        return Tg_ * Ar_moleculardynamics::YPSILON / Ar_moleculardynamics::KB;
    }

    void Ar_moleculardynamics::Move_Atoms()
    {
        // calculate temperture
        auto Uk = 0.0;

        for (auto n = 0; n < NumAtom; n++) {
            auto const v2 = norm2(VX[n], VY[n], VZ[n]);
            Uk += v2;
        }

        // 運動エネルギーの計算
        Uk *= 0.5;

        // 温度の計算
        Tc_ = Uk / (1.5 * static_cast<double>(NumAtom));

        auto const s = std::sqrt((Tg_ + Ar_moleculardynamics::ALPHA * (Tc_ - Tg_)) / Tc_);
		
        switch (MD_iter_) {
        case 1:
            // update the coordinates by the second order Euler method
            // 最初のステップだけ修正Euler法で時間発展
			tbb::parallel_for(
				0,
				NumAtom,
				1,
				[this, s](std::int32_t n) {
				X1[n] = X_[n];
				Y1[n] = Y_[n];
				Z1[n] = Z_[n];

				// scaling of velocity
				VX[n] *= s;
				VY[n] *= s;
				VZ[n] *= s;

				// update coordinates and velocity
				X_[n] += Ar_moleculardynamics::DT * VX[n] + 0.5 * FX[n] * dt2;
				Y_[n] += Ar_moleculardynamics::DT * VY[n] + 0.5 * FY[n] * dt2;
				Z_[n] += Ar_moleculardynamics::DT * VZ[n] + 0.5 * FZ[n] * dt2;

				VX[n] += Ar_moleculardynamics::DT * FX[n];
				VY[n] += Ar_moleculardynamics::DT * FY[n];
				VZ[n] += Ar_moleculardynamics::DT * FZ[n];
			},
				tbb::auto_partitioner());
            break;

        default:
            // update the coordinates by the Verlet method

			tbb::parallel_for(
				0,
				NumAtom,
				1,
				[this, s](std::int32_t n) {
				std::array<double, 3> rtmp;
				rtmp[0] = X_[n];
				rtmp[1] = Y_[n];
				rtmp[2] = Z_[n];

				// update coordinates and velocity
				// Verlet法の座標更新式において速度成分を抜き出し、その部分をスケールする
//#ifdef NVT
				X_[n] += s * (X_[n] - X1[n]) + FX[n] * dt2;
				Y_[n] += s * (Y_[n] - Y1[n]) + FY[n] * dt2;
				Z_[n] += s * (Z_[n] - Z1[n]) + FZ[n] * dt2;
//#elif NVE
				//X_[n] = 2.0 * X_[n] - X1[n] + FX[n] * dt2;
				//Y_[n] = 2.0 * Y_[n] - Y1[n] + FY[n] * dt2;
				//Z_[n] = 2.0 * Z_[n] - Z1[n] + FZ[n] * dt2;
//#endif
				VX[n] = 0.5 * (X_[n] - X1[n]) / Ar_moleculardynamics::DT;
				VY[n] = 0.5 * (Y_[n] - Y1[n]) / Ar_moleculardynamics::DT;
				VZ[n] = 0.5 * (Z_[n] - Z1[n]) / Ar_moleculardynamics::DT;

				X1[n] = rtmp[0];
				Y1[n] = rtmp[1];
				Z1[n] = rtmp[2];
			},
				tbb::auto_partitioner());
            break;
        }

        // consider the periodic boundary condination
        // セルの外側に出たら座標をセル内に戻す
        tbb::parallel_for(
            0,
            NumAtom,
            1,
            [this](std::int32_t n) {
			if (X_[n] > periodiclen) {
                X_[n] -= periodiclen;
                X1[n] -= periodiclen;
            }
            else if (X_[n] < 0.0) {
                X_[n] += periodiclen;
                X1[n] += periodiclen;
            }
            if (Y_[n] > periodiclen) {
                Y_[n] -= periodiclen;
                Y1[n] -= periodiclen;
            }
            else if (Y_[n] < 0.0) {
                Y_[n] += periodiclen;
                Y1[n] += periodiclen;
            }
            if (Z_[n] > periodiclen) {
                Z_[n] -= periodiclen;
                Z1[n] -= periodiclen;
            }
            else if (Z_[n] < 0.0) {
                Z_[n] += periodiclen;
                Z1[n] += periodiclen;
            }
        },
            tbb::auto_partitioner());

        // 繰り返し回数と時間を増加
        t = static_cast<double>(MD_iter_) * Ar_moleculardynamics::DT;
        MD_iter_++;
    }

    void Ar_moleculardynamics::recalc()
    {
		t = 0.0;
		MD_iter_ = 1;

        MD_initPos();
        MD_initVel();
    }

	void Ar_moleculardynamics::setScale(double scale)
	{
		scale_ = scale;
		lat_ = std::pow(2.0, 2.0 / 3.0) * scale_;

		recalc();

		periodiclen = lat_ * static_cast<double>(Nc);
	}

	void Ar_moleculardynamics::setTgiven(double Tgiven)
	{
		Tg_ = Tgiven * Ar_moleculardynamics::KB / Ar_moleculardynamics::YPSILON;
	}

    // #endregion publicメンバ関数

	// #region privateメンバ関数

	void Ar_moleculardynamics::MD_initPos()
	{
		double sx, sy, sz;
		auto n = 0;

		for (auto i = 0; i < Nc; i++) {
			for (auto j = 0; j < Nc; j++) {
				for (auto k = 0; k < Nc; k++) {
				    // 基本セルをコピーする
				    sx = static_cast<double>(i) * lat_;
				    sy = static_cast<double>(j) * lat_;
				    sz = static_cast<double>(k) * lat_;

				    // 基本セル内には4つの原子がある
				    X_[n] = sx;
				    Y_[n] = sy;
				    Z_[n] = sz;
				    n++;

				    X_[n] = 0.5 * lat_ + sx;
				    Y_[n] = 0.5 * lat_ + sy;
				    Z_[n] = sz;
				    n++;

				    X_[n] = sx;
				    Y_[n] = 0.5 * lat_ + sy;
				    Z_[n] = 0.5 * lat_ + sz;
				    n++;

				    X_[n] = 0.5 * lat_ + sx;
				    Y_[n] = sy;
				    Z_[n] = 0.5 * lat_ + sz;
				    n++;
				}
            }
        }

		NumAtom = n;

		// move the center of mass to the origin
		// 系の重心を座標系の原点とする
		sx = 0.0;
		sy = 0.0;
		sz = 0.0;

		for (auto n = 0; n < NumAtom; n++) {
			sx += X_[n];
			sy += Y_[n];
			sz += Z_[n];
		}

		sx /= static_cast<double>(NumAtom);
		sy /= static_cast<double>(NumAtom);
		sz /= static_cast<double>(NumAtom);

		for (auto n = 0; n < NumAtom; n++) {
			X_[n] -= sx;
			Y_[n] -= sy;
			Z_[n] -= sz;
		}
	}

	void Ar_moleculardynamics::MD_initVel()
	{
		auto const v = std::sqrt(3.0 * Tg_);

		myrandom::MyRand mr(-1.0, 1.0);

		for (auto n = 0; n < NumAtom; n++) {
			double rndX = mr.myrand();
			double rndY = mr.myrand();
			double rndZ = mr.myrand();
			double tmp = 1.0 / std::sqrt(norm2(rndX, rndY, rndZ));
			rndX *= tmp;
			rndY *= tmp;
			rndZ *= tmp;

			// 方向はランダムに与える
			VX[n] = v * rndX;
			VY[n] = v * rndY;
			VZ[n] = v * rndZ;
		}

		auto sx = 0.0;
		auto sy = 0.0;
		auto sz = 0.0;

		for (auto n = 0; n < NumAtom; n++) {
			sx += VX[n];
			sy += VY[n];
			sz += VZ[n];
		}

		sx /= static_cast<double>(NumAtom);
		sy /= static_cast<double>(NumAtom);
		sz /= static_cast<double>(NumAtom);

		// 重心の並進運動を避けるために、速度の和がゼロになるように補正
		for (auto n = 0; n < NumAtom; n++) {
			VX[n] -= sx;
			VY[n] -= sy;
			VZ[n] -= sz;
		}
	}

	double Ar_moleculardynamics::norm2(double x, double y, double z) const
	{
		return (x * x + y * y + z * z);
	}

	// #endregion privateメンバ関数
}
