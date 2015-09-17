/*! \file Ar_moleculardynamics.h
	\brief アルゴンに対して、分子動力学シミュレーションを行うクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#include "DXUT.h"
#include "Ar_moleculardynamics.h"
#include <tbb/parallel_for.h>							// for tbb::parallel_for
#include <tbb/partitioner.h>							// for tbb::auto_partitioner
#include <boost/simd/include/functions/divides.hpp>
#include <boost/simd/include/functions/plus.hpp>
#include <boost/simd/include/functions/sum.hpp>
#include <boost/simd/include/functions/multiplies.hpp>

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
		C([this] { return std::cref(C_); }, nullptr),
		MD_iter([this] { return MD_iter_; }, nullptr),
		Nc([this] { return Nc_; }, nullptr),
		NumAtom([this] { return NumAtom_; }, nullptr),
		periodiclen([this] { return periodiclen_; }, nullptr),
		useavx([this] { return useavx_; }, nullptr),
		X([this] { return std::cref(X_); }, nullptr),
		Y([this] { return std::cref(Y_); }, nullptr),
		Z([this] { return std::cref(Z_); }, nullptr),
		C_(Nc_ * Nc_ * Nc_ * 4),
		C1_(Nc_ * Nc_ * Nc_ * 4),
		dt2_(DT * DT),
		F_(Nc_ * Nc_ * Nc_ * 4),
		FX_(Nc_ * Nc_ * Nc_ * 4),
		FY_(Nc_ * Nc_ * Nc_ * 4),
		FZ_(Nc_ * Nc_ * Nc_ * 4),
		rc2_(rc_ * rc_),
		rcm6_(std::pow(rc_, -6.0)),
		rcm12_(std::pow(rc_, -12.0)),
		Tg_(Ar_moleculardynamics::FIRSTTEMP * Ar_moleculardynamics::KB / Ar_moleculardynamics::YPSILON),
		useavx_(true),
		Vrc_(4.0 * (rcm12_ - rcm6_)),
		V_(Nc_ * Nc_ * Nc_ * 4),
		VX_(Nc_ * Nc_ * Nc_ * 4),
		VY_(Nc_ * Nc_ * Nc_ * 4),
		VZ_(Nc_ * Nc_ * Nc_ * 4),
		X_(Nc_ * Nc_ * Nc_ * 4),
		X1_(Nc_ * Nc_ * Nc_ * 4),
		Y_(Nc_ * Nc_ * Nc_ * 4),
		Y1_(Nc_ * Nc_ * Nc_ * 4),
		Z_(Nc_ * Nc_ * Nc_ * 4),
		Z1_(Nc_ * Nc_ * Nc_ * 4)
	{
		// initalize parameters
		lat_ = std::pow(2.0, 2.0 / 3.0) * scale_;

        recalc();

		periodiclen_ = lat_ * static_cast<double>(Nc_);
	}

	// #endregion コンストラクタ

    // #region publicメンバ関数

	template <>
	void Ar_moleculardynamics::Calc_Forces<UseAVX::True>()
	{
		for (auto n = 0; n < NumAtom_; n++) {
			F_[n] = Ar_moleculardynamics::pack_t(0.0, 0.0, 0.0, 0.0);
		}

		tbb::parallel_for(
				0,
				NumAtom_,
				1,
				[this](std::int32_t n) {
			for (auto m = 0; m < NumAtom_; m++) {

				// ±ncp分のセル内の原子との相互作用を計算
				for (auto i = -ncp; i <= ncp; i++) {
					for (auto j = -ncp; j <= ncp; j++) {
						for (auto k = -ncp; k <= ncp; k++) {
							auto const sx = static_cast<double>(i)* periodiclen_;
							auto const sy = static_cast<double>(j)* periodiclen_;
							auto const sz = static_cast<double>(k)* periodiclen_;

							// 自分自身との相互作用を排除
							if (n != m || i != 0 || j != 0 || k != 0) {
								auto dx = C_[n][0] - (C_[m][0] + sx);
								auto dy = C_[n][1] - (C_[m][1] + sy);
								auto dz = C_[n][2] - (C_[m][2] + sz);

								auto const r2 = norm2(dx, dy, dz);
								// 打ち切り距離内であれば計算
								if (r2 <= rc2_) {
									auto const r = std::sqrt(r2);
									auto const rm6 = 1.0 / (r2 * r2 * r2);
									auto const rm7 = rm6 / r;
									auto const rm12 = rm6 * rm6;
									auto const rm13 = rm12 / r;

									auto const Fr = 48.0 * rm13 - 24.0 * rm7;
									auto const Frdivr = Fr / r;

									Ar_moleculardynamics::pack_t fr(Frdivr, Frdivr, Frdivr, 0.0);

									Ar_moleculardynamics::pack_t p(dx, dy, dz, 0.0);

									F_[n] += p * fr;
								}
							}
						}
					}
				}
			}
		},
		tbb::auto_partitioner());
	}

	template <>
	void Ar_moleculardynamics::Calc_Forces<UseAVX::False>()
    {
        for (auto n = 0; n < NumAtom_; n++) {
            FX_[n] = 0.0;
            FY_[n] = 0.0;
            FZ_[n] = 0.0;
        }

		tbb::parallel_for(
			0,
			NumAtom_,
			1,
			[this](std::int32_t n) {
			for (auto m = 0; m < NumAtom_; m++) {

				// ±ncp分のセル内の原子との相互作用を計算
				for (auto i = -ncp; i <= ncp; i++) {
					for (auto j = -ncp; j <= ncp; j++) {
						for (auto k = -ncp; k <= ncp; k++) {
							auto const sx = static_cast<double>(i) * periodiclen_;
							auto const sy = static_cast<double>(j) * periodiclen_;
							auto const sz = static_cast<double>(k) * periodiclen_;

							// 自分自身との相互作用を排除
							if (n != m || i != 0 || j != 0 || k != 0) {
								auto const dx = X_[n] - (X_[m] + sx);
								auto const dy = Y_[n] - (Y_[m] + sy);
								auto const dz = Z_[n] - (Z_[m] + sz);

								auto const r2 = norm2(dx, dy, dz);
								// 打ち切り距離内であれば計算
								if (r2 <= rc2_) {
									auto const r = std::sqrt(r2);
									auto const rm6 = 1.0 / (r2 * r2 * r2);
									auto const rm7 = rm6 / r;
									auto const rm12 = rm6 * rm6;
									auto const rm13 = rm12 / r;

									auto const Fr = 48.0 * rm13 - 24.0 * rm7;
									auto const Frdivr = Fr / r;

									FX_[n] += dx * Frdivr;
									FY_[n] += dy * Frdivr;
									FZ_[n] += dz * Frdivr;
								}
							}
						}
					}
				}
			}
		},
		tbb::auto_partitioner());
	}
	
	double Ar_moleculardynamics::getDeltat() const
	{
		return Ar_moleculardynamics::TAU * t_ * 1.0E+12;
	}

    float Ar_moleculardynamics::getForce(std::int32_t n) const
    {
		auto const force =
			useavx_ ? std::sqrt(norm2(F_[n][0], F_[n][1], F_[n][2])) : std::sqrt(norm2(FX_[n], FY_[n], FZ_[n]));

		return static_cast<float>(force);
    }

	double Ar_moleculardynamics::getLatticeconst() const
	{
		return Ar_moleculardynamics::SIGMA * lat_ * 1.0E+10;
	}

	double Ar_moleculardynamics::getPeriodiclen() const
	{
		return Ar_moleculardynamics::SIGMA * periodiclen_ * 1.0E+10;
	}

    double Ar_moleculardynamics::getTcalc() const
    {
		return Ar_moleculardynamics::YPSILON / Ar_moleculardynamics::KB * Tc_;
    }

    double Ar_moleculardynamics::getTgiven() const
    {
		return Ar_moleculardynamics::YPSILON / Ar_moleculardynamics::KB * Tg_;
    }

	template<>
	void Ar_moleculardynamics::Move_Atoms<UseAVX::True>()
	{
		// calculate temperture
		auto Uk = 0.0;

		for (auto n = 0; n < NumAtom_; n++) {
			auto const v2 = norm2(V_[n][0], V_[n][1], V_[n][2]);
			Uk += v2;
		}

		// 運動エネルギーの計算
		Uk *= 0.5;

		// 温度の計算
		Tc_ = Uk / (1.5 * static_cast<double>(NumAtom_));

		auto const s = std::sqrt((Tg_ + Ar_moleculardynamics::ALPHA * (Tc_ - Tg_)) / Tc_);

		switch (MD_iter_) {
		case 1:
			// update the coordinates by the second order Euler method
			// 最初のステップだけ修正Euler法で時間発展
			tbb::parallel_for(
				0,
				NumAtom_,
				1,
				[this, s](std::int32_t n) {
				C1_[n] = C_[n];

				// scaling of velocity
				V_[n] *= s;

				// update coordinates and velocity
				C_[n] += Ar_moleculardynamics::DT * V_[n] + 0.5 * F_[n] * dt2_;

				V_[n] += Ar_moleculardynamics::DT * F_[n];

			},
				tbb::auto_partitioner());
			break;

		default:
			// update the coordinates by the Verlet method

			tbb::parallel_for(
				0,
				NumAtom_,
				1,
				[this, s](std::int32_t n) {
				auto const rtmp = C_[n];

				// update coordinates and velocity
				// Verlet法の座標更新式において速度成分を抜き出し、その部分をスケールする
				C_[n] += s * (C_[n] - C1_[n]) + F_[n] * dt2_;

				V_[n] = 0.5 * (C_[n] - C1_[n]) / Ar_moleculardynamics::DT;

				C1_[n] = rtmp;
			},
				tbb::auto_partitioner());
			break;
		}

		// consider the periodic boundary condination
		// セルの外側に出たら座標をセル内に戻す
		tbb::parallel_for(
			0,
			NumAtom_,
			1,
			[this](std::int32_t n) {
			if (C_[n][0] > periodiclen_) {
				C_[n][0] -= periodiclen_;
				C1_[n][0] -= periodiclen_;
			}
			else if (C_[n][0] < 0.0) {
				C_[n][0] += periodiclen_;
				C1_[n][0] += periodiclen_;
			}
			if (C_[n][1] > periodiclen_) {
				C_[n][1] -= periodiclen_;
				C1_[n][1] -= periodiclen_;
			}
			else if (C_[n][1] < 0.0) {
				C_[n][1] += periodiclen_;
				C1_[n][1] += periodiclen_;
			}
			if (C_[n][2] > periodiclen_) {
				C_[n][2] -= periodiclen_;
				C1_[n][2] -= periodiclen_;
			}
			else if (C_[n][2] < 0.0) {
				C_[n][2] += periodiclen_;
				C1_[n][2] += periodiclen_;
			}
		},
			tbb::auto_partitioner());
		
		// 繰り返し回数と時間を増加
		t_ = static_cast<double>(MD_iter_)* Ar_moleculardynamics::DT;
		MD_iter_++;
	}

	template<>
	void Ar_moleculardynamics::Move_Atoms<UseAVX::False>()
    {
        // calculate temperture
        auto Uk = 0.0;

        for (auto n = 0; n < NumAtom_; n++) {
            auto const v2 = norm2(VX_[n], VY_[n], VZ_[n]);
            Uk += v2;
        }

        // 運動エネルギーの計算
        Uk *= 0.5;

        // 温度の計算
        Tc_ = Uk / (1.5 * static_cast<double>(NumAtom_));

        auto const s = std::sqrt((Tg_ + Ar_moleculardynamics::ALPHA * (Tc_ - Tg_)) / Tc_);
		
        switch (MD_iter_) {
        case 1:
            // update the coordinates by the second order Euler method
            // 最初のステップだけ修正Euler法で時間発展
			tbb::parallel_for(
				0,
				NumAtom_,
				1,
				[this, s](std::int32_t n) {
				X1_[n] = X_[n];
				Y1_[n] = Y_[n];
				Z1_[n] = Z_[n];

				// scaling of velocity
				VX_[n] *= s;
				VY_[n] *= s;
				VZ_[n] *= s;

				// update coordinates and velocity
				X_[n] += Ar_moleculardynamics::DT * VX_[n] + 0.5 * FX_[n] * dt2_;
				Y_[n] += Ar_moleculardynamics::DT * VY_[n] + 0.5 * FY_[n] * dt2_;
				Z_[n] += Ar_moleculardynamics::DT * VZ_[n] + 0.5 * FZ_[n] * dt2_;

				VX_[n] += Ar_moleculardynamics::DT * FX_[n];
				VY_[n] += Ar_moleculardynamics::DT * FY_[n];
				VZ_[n] += Ar_moleculardynamics::DT * FZ_[n];
			},
				tbb::auto_partitioner());
            break;

        default:
            // update the coordinates by the Verlet method

			tbb::parallel_for(
				0,
				NumAtom_,
				1,
				[this, s](std::int32_t n) {
				std::array<double, 3> rtmp;
				rtmp[0] = X_[n];
				rtmp[1] = Y_[n];
				rtmp[2] = Z_[n];

				// update coordinates and velocity
				// Verlet法の座標更新式において速度成分を抜き出し、その部分をスケールする
				X_[n] += s * (X_[n] - X1_[n]) + FX_[n] * dt2_;
				Y_[n] += s * (Y_[n] - Y1_[n]) + FY_[n] * dt2_;
				Z_[n] += s * (Z_[n] - Z1_[n]) + FZ_[n] * dt2_;

				VX_[n] = 0.5 * (X_[n] - X1_[n]) / Ar_moleculardynamics::DT;
				VY_[n] = 0.5 * (Y_[n] - Y1_[n]) / Ar_moleculardynamics::DT;
				VZ_[n] = 0.5 * (Z_[n] - Z1_[n]) / Ar_moleculardynamics::DT;

				X1_[n] = rtmp[0];
				Y1_[n] = rtmp[1];
				Z1_[n] = rtmp[2];
			},
				tbb::auto_partitioner());
            break;
        }

        // consider the periodic boundary condination
        // セルの外側に出たら座標をセル内に戻す
        tbb::parallel_for(
            0,
            NumAtom_,
            1,
            [this](std::int32_t n) {
			if (X_[n] > periodiclen_) {
                X_[n] -= periodiclen_;
                X1_[n] -= periodiclen_;
            }
            else if (X_[n] < 0.0) {
                X_[n] += periodiclen_;
                X1_[n] += periodiclen_;
            }
            if (Y_[n] > periodiclen_) {
                Y_[n] -= periodiclen_;
                Y1_[n] -= periodiclen_;
            }
            else if (Y_[n] < 0.0) {
                Y_[n] += periodiclen_;
                Y1_[n] += periodiclen_;
            }
            if (Z_[n] > periodiclen_) {
                Z_[n] -= periodiclen_;
                Z1_[n] -= periodiclen_;
            }
            else if (Z_[n] < 0.0) {
                Z_[n] += periodiclen_;
                Z1_[n] += periodiclen_;
            }
        },
            tbb::auto_partitioner());

        // 繰り返し回数と時間を増加
        t_ = static_cast<double>(MD_iter_) * Ar_moleculardynamics::DT;
        MD_iter_++;
    }

    void Ar_moleculardynamics::recalc()
    {
		t_ = 0.0;
		MD_iter_ = 1;

		if (useavx_) {
			MD_initPos<UseAVX::True>();
			MD_initVel<UseAVX::True>();
		}
		else {
			MD_initPos<UseAVX::False>();
			MD_initVel<UseAVX::False>();
		}		
    }

	void Ar_moleculardynamics::setNc(std::int32_t Nc)
	{
		Nc_ = Nc;
		F_.resize(Nc_ * Nc_ * Nc_ * 4);
		
		FX_.resize(Nc_ * Nc_ * Nc_ * 4);
		FY_.resize(Nc_ * Nc_ * Nc_ * 4);
		FZ_.resize(Nc_ * Nc_ * Nc_ * 4);

		V_.resize(Nc_ * Nc_ * Nc_ * 4);

		VX_.resize(Nc_ * Nc_ * Nc_ * 4);
		VY_.resize(Nc_ * Nc_ * Nc_ * 4);
		VZ_.resize(Nc_ * Nc_ * Nc_ * 4);

		C_.resize(Nc_ * Nc_ * Nc_ * 4);
		C1_.resize(Nc_ * Nc_ * Nc_ * 4);

		X_.resize(Nc_ * Nc_ * Nc_ * 4);
		X1_.resize(Nc_ * Nc_ * Nc_ * 4);
		Y_.resize(Nc_ * Nc_ * Nc_ * 4);
		Y1_.resize(Nc_ * Nc_ * Nc_ * 4);
		Z_.resize(Nc_ * Nc_ * Nc_ * 4);
		Z1_.resize(Nc_ * Nc_ * Nc_ * 4);

		ModLattice();
	}

	void Ar_moleculardynamics::setScale(double scale)
	{
		scale_ = scale;
		ModLattice();
	}

	void Ar_moleculardynamics::setTgiven(double Tgiven)
	{
		Tg_ = Tgiven * Ar_moleculardynamics::KB / Ar_moleculardynamics::YPSILON;
	}

    // #endregion publicメンバ関数

	// #region privateメンバ関数

    bool Ar_moleculardynamics::availableAVX() const
    {
#if (_MSC_FULL_VER >= 160040219)
        std::array<std::int32_t, 4> cpuInfo = { 0 };
        ::__cpuid(cpuInfo.data(), 1);

        auto const osUsesXSAVE_XRSTORE = cpuInfo[2] & (1 << 27) || false;
        auto const cpuAVXSuport = cpuInfo[2] & (1 << 28) || false;

        if (osUsesXSAVE_XRSTORE && cpuAVXSuport)
        {
            // Check if the OS will save the YMM registers
            auto const xcrFeatureMask = ::_xgetbv(_XCR_XFEATURE_ENABLED_MASK);
            return (xcrFeatureMask & 0x6) || false;
        }
#endif
        return false;
    }

	void Ar_moleculardynamics::ModLattice()
	{
		lat_ = std::pow(2.0, 2.0 / 3.0) * scale_;
		recalc();
		periodiclen_ = lat_ * static_cast<double>(Nc_);
	}

	double Ar_moleculardynamics::norm2(double x, double y, double z) const
	{
		return (x * x + y * y + z * z);
	}

	// #endregion privateメンバ関数
	
	// #region templateメンバ関数の実体化

	template void Ar_moleculardynamics::Calc_Forces<UseAVX::True>();
	template void Ar_moleculardynamics::Calc_Forces<UseAVX::False>();
	template void Ar_moleculardynamics::Move_Atoms<UseAVX::True>();
	template void Ar_moleculardynamics::Move_Atoms<UseAVX::False>();

	// #endregion templateメンバ関数の実体化
}
