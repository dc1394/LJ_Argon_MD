/*! \file Ar_moleculardynamics.h
	\brief アルゴンに対して、分子動力学シミュレーションを行うクラスの宣言

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#ifndef _AR_MOLECULARDYNAMICS_H_
#define _AR_MOLECULARDYNAMICS_H_

#pragma once

#include "../myrandom/myrand.h"
#include "../utility/property.h"
#include "useavx.h"
#include <array>							// for std::array
#include <cmath>							// for std::sqrt
#include <cstdint>							// for std::int32_t
#include <vector>							// for std::vector
#include <dvec.h>
#include <boost/simd/memory/allocator.hpp>  // for boost::simd::allocator
#include <boost/simd/sdk/simd/pack.hpp>		// for boost::simd::pack

namespace moleculardynamics {
	using namespace utility;

    //! A class.
    /*!
        アルゴンに対して、分子動力学シミュレーションを行うクラス
    */
    class Ar_moleculardynamics final {
		// #region 型エイリアス

		using pack_t = boost::simd::pack<double>;
		
		// #endregion 型エイリアス

        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            コンストラクタ
        */
        Ar_moleculardynamics();

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~Ar_moleculardynamics() = default;

        // #endregion コンストラクタ・デストラクタ

		// #region publicメンバ関数

        template <UseAVX U>
        //! A public member function.
        /*!
            原子に働く力を計算する
        */
        void Calc_Forces();
		
		//! A public member function (constant).
		/*!
			シミュレーションを開始してからの経過時間を求める
		*/
		double getDeltat() const;

		//! A public member function (constant).
        /*!
            n番目の原子に働く力を求める
        */
        float getForce(std::int32_t n) const;

		//! A public member function (constant).
		/*!
			格子定数を求める
		*/
		double getLatticeconst() const;
		
		//! A public member function (constant).
		/*!
			周期境界条件の長さを求める
		*/
		double getPeriodiclen() const;

        //! A public member function (constant).
        /*!
            計算された温度の絶対温度を求める
        */
        double getTcalc() const;

        //! A public member function (constant).
        /*!
            与えた温度の絶対温度を求める
        */
        double getTgiven() const;
				
		template <UseAVX U>
        //! A public member function.
        /*!
            原子を移動させる
        */
        void Move_Atoms();
		
		//! A oublic member function.
		/*!
			再計算する
		*/
		void recalc();

		//! A public member function.
		/*!
			スーパーセルの大きさを設定する
			\param Nc_ スーパーセルの大きさ
		*/
		void setNc(std::int32_t Nc);

		//! A public member function.
		/*!
			格子定数のスケールを設定する
			\param scale 設定する格子定数のスケール
		*/
		void setScale(double scale);

		//! A public member function.
		/*!
			温度を設定する
			\param Tgiven 設定する温度（絶対温度）
		*/
		void setTgiven(double Tgiven);

        // #endregion publicメンバ関数

		// #region privateメンバ関数

	private:
        //! A private member function (constant).
        /*!
            AVX命令が使えるかどうか
        */
        bool availableAVX() const;

		template <UseAVX U>
		//! A private member function.
		/*!
			原子の初期位置を決める
		*/
		void MD_initPos();

		template <UseAVX U>
		//! A private member function.
		/*!
			原子の初期速度を決める
		*/
		void MD_initVel();

		//! A private member function.
		/*!
			格子定数が変更されたときに呼ばれる
		*/
		void ModLattice();

		//! A private member function (constant).
		/*!
			ノルムの二乗を求める
			\param x x座標
			\param y y座標
			\param z z座標
			\return ノルムの二乗
		*/
		double norm2(double x, double y, double z) const;

		// #endregion privateメンバ関数

		// #region プロパティ

	public:
		//! A property.
		/*!
			n番目の原子の座標へのプロパティ
		*/
		Property<std::vector<double, boost::simd::allocator<double>> const &> const C;

		//! A property.
		/*!
			MDのステップ数へのプロパティ
		*/
		Property<std::int32_t> const MD_iter;
		
		//! A property.
		/*!
			スーパーセルの個数へのプロパティ
		*/
		Property<std::int32_t> const Nc;

		//! A property.
		/*!
			原子数へのプロパティ
		*/
		Property<std::int32_t> const NumAtom;

		//! A property.
		/*!
			格子定数へのプロパティ
		*/
		Property<double> const periodiclen;

		//! A property.
		/*!
			AVX命令が使えるかどうかへのプロパティ
		*/
		Property<bool> const useavx;

		//! A property.
		/*!
			n番目の原子のx座標へのプロパティ
		*/
		Property<std::vector<double> const &> const X;

		//! A property.
		/*!
			n番目の原子のy座標へのプロパティ
		*/
		Property<std::vector<double> const &> const Y;

		//! A property.
		/*!
			n番目の原子のz座標へのプロパティ
		*/
		Property<std::vector<double> const &> const Z;

		// #endregion プロパティ

        // #region メンバ変数

	public:
		//! A private member variable (constant).
		/*!
			初期のスーパーセルの個数
		*/
		static auto const FIRSTNC = 4;

		//! A private member variable (constant).
		/*!
			初期の格子定数のスケール
		*/
		static double const FIRSTSCALE;

		//! A private member variable (constant).
		/*!
			初期温度（絶対温度）
		*/
		static double const FIRSTTEMP;

	private:
		//! A private member variable (constant).
		/*!
			Woodcockの温度スケーリングの係数
		*/
		static double const ALPHA;

		//! A private member variable (constant).
		/*!
			アボガドロ定数
		*/
		static double const AVOGADRO_CONSTANT;

		//! A private member variable (constant).
		/*!
			時間刻みΔt
		*/
		static double const DT;
				
		//! A private member variable (constant).
        /*!
            ボルツマン定数
        */
        static double const KB;
        
        //! A private member variable (constant).
        /*!
            アルゴン原子に対するσ
        */
        static double const SIGMA;

		//! A private member variable (constant).
		/*!
			アルゴン原子に対するτ
		*/
		static double const TAU;

        //! A private member variable (constant).
        /*!
            アルゴン原子に対するε
        */
        static double const YPSILON;
		
		//! A private member variable (constant).
		/*!
			スーパーセルの個数
		*/
		std::int32_t Nc_ = Ar_moleculardynamics::FIRSTNC;

		//! A private member variable.
		/*!
			n個目の原子の座標
		*/
		std::vector<double, boost::simd::allocator<double>> C_;

		//! A private member variable.
		/*!
			n個目の原子の前の計算の座標
		*/
		std::vector<double, boost::simd::allocator<double>> C1_;
		
		//! A private member variable (constant).
		/*!
			時間刻みの二乗
		*/
		double const dt2_;

		//! A private member variable.
		/*!
			n個目の原子に働く力のx成分
		*/
		std::vector<double, boost::simd::allocator<double>> F_;

		//! A private member variable.
		/*!
			n個目の原子に働く力のx成分
		*/
		std::vector<double> FX_;

		//! A private member variable.
		/*!
			n個目の原子に働く力のy成分
		*/
		std::vector<double> FY_;

		//! A private member variable.
		/*!
			n個目の原子に働く力のz成分
		*/
		std::vector<double> FZ_;
		
		//! A private member variable.
		/*!
			格子定数
		*/
		double lat_;
		
		//! A private member variable.
		/*!
			MDのステップ数
		*/
		std::int32_t MD_iter_;

		//! A private member variable (constant).
		/*!
			相互作用を計算するセルの個数
		*/
		std::int32_t const ncp = 2;
		
		//! A private member variable.
		/*!
			原子数
		*/
		std::int32_t NumAtom_;
		
		//! A private member variable.
		/*!
			周期境界条件の長さ
		*/
		double periodiclen_;

		//! A private member variable (constant).
		/*!
			カットオフ半径
		*/
		double const rc_ = 2.5;

		//! A private member variable (constant).
		/*!
			カットオフ半径の2乗
		*/
		double const rc2_;

		//! A private member variable (constant).
		/*!
			カットオフ半径の逆数の6乗
		*/
		double const rcm6_;

		//! A private member variable (constant).
		/*!
			カットオフ半径の逆数の12乗
		*/
		double const rcm12_;

		//! A private member variable.
		/*!
			格子定数のスケーリングの定数
		*/
		double scale_ = Ar_moleculardynamics::FIRSTSCALE;
		
		//! A private member variable.
		/*!
			時間	
		*/
		double t_;

		//! A private member variable.
		/*!
			計算された温度Tcalc
		*/
        double Tc_;

		//! A private member variable.
		/*!
			与える温度Tgiven
		*/
        double Tg_;

		//! A private member variable (constant).
		/*!
			AVX命令を使うかどうか
		*/
		bool const useavx_;

		//! A private member variable (constant).
		/*!
			ポテンシャルエネルギーの打ち切り
		*/
        double const Vrc_;

		//! A private member variable.
		/*!
			n個目の原子の速度
		*/
		std::vector<double, boost::simd::allocator<double>> V_;

		//! A private member variable.
		/*!
			n個目の原子の速度のx成分
		*/
		std::vector<double> VX_;

		//! A private member variable.
		/*!
			n個目の原子の速度のy成分
		*/
		std::vector<double> VY_;

		//! A private member variable.
		/*!
			n個目の原子の速度のz成分
		*/
		std::vector<double> VZ_;

		//! A private member variable.
		/*!
			n個目の原子のx座標
		*/
        std::vector<double> X_;

		//! A private member variable.
		/*!
			n個目の原子の前の計算のx座標
		*/
		std::vector<double> X1_;
		
		//! A private member variable.
		/*!
			n個目の原子のy座標
		*/
		std::vector<double> Y_;

		//! A private member variable.
		/*!
			n個目の原子の前の計算のy座標
		*/
		std::vector<double> Y1_;
		
		//! A private member variable.
		/*!
			n個目の原子のz座標
		*/
		std::vector<double> Z_;
        
		//! A private member variable.
		/*!
			n個目の原子の前の計算のz座標
		*/
		std::vector<double> Z1_;

		// #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        Ar_moleculardynamics(Ar_moleculardynamics const &) = delete;

        //! A private member function (deleted).
        /*!
			operator=()の宣言（禁止）
			\param コピー元のオブジェクト（未使用）
			\return コピー元のオブジェクト
        */
        Ar_moleculardynamics & operator=(Ar_moleculardynamics const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };

	// #region privateメンバ関数
	
	template<>
	inline void Ar_moleculardynamics::MD_initPos<UseAVX::True>()
	{
		double sx, sy, sz;
		auto n = 0;

		for (auto i = 0; i < Nc_; i++) {
			for (auto j = 0; j < Nc_; j++) {
				for (auto k = 0; k < Nc_; k++) {
					// 基本セルをコピーする
					sx = static_cast<double>(i) * lat_;
					sy = static_cast<double>(j) * lat_;
					sz = static_cast<double>(k) * lat_;

					// 基本セル内には4つの原子がある
					C_[n][0] = sx;
					C_[n][1] = sy;
					C_[n][2] = sz;
					n++;

					C_[n][0] = 0.5 * lat_ + sx;
					C_[n][1] = 0.5 * lat_ + sy;
					C_[n][2] = sz;
					n++;

					C_[n][0] = sx;
					C_[n][1] = 0.5 * lat_ + sy;
					C_[n][2] = 0.5 * lat_ + sz;
					n++;

					C_[n][0] = 0.5 * lat_ + sx;
					C_[n][1] = sy;
					C_[n][2] = 0.5 * lat_ + sz;
					n++;
				}
			}
		}

		NumAtom_ = n;

		// move the center of mass to the origin
		// 系の重心を座標系の原点とする
		sx = 0.0;
		sy = 0.0;
		sz = 0.0;

		for (auto n = 0; n < NumAtom_; n++) {
			sx += C_[n][0];
			sy += C_[n][1];
			sz += C_[n][2];
		}

		sx /= static_cast<double>(NumAtom_);
		sy /= static_cast<double>(NumAtom_);
		sz /= static_cast<double>(NumAtom_);

		for (auto n = 0; n < NumAtom_; n++) {
			C_[n][0] -= sx;
			C_[n][1] -= sy;
			C_[n][2] -= sz;
		}
	}

	template<>
	inline void Ar_moleculardynamics::MD_initPos<UseAVX::False>()
	{
		double sx, sy, sz;
		auto n = 0;

		for (auto i = 0; i < Nc_; i++) {
			for (auto j = 0; j < Nc_; j++) {
				for (auto k = 0; k < Nc_; k++) {
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

		NumAtom_ = n;

		// move the center of mass to the origin
		// 系の重心を座標系の原点とする
		sx = 0.0;
		sy = 0.0;
		sz = 0.0;

		for (auto n = 0; n < NumAtom_; n++) {
			sx += X_[n];
			sy += Y_[n];
			sz += Z_[n];
		}

		sx /= static_cast<double>(NumAtom_);
		sy /= static_cast<double>(NumAtom_);
		sz /= static_cast<double>(NumAtom_);

		for (auto n = 0; n < NumAtom_; n++) {
			X_[n] -= sx;
			Y_[n] -= sy;
			Z_[n] -= sz;
		}
	}

	template<>
	inline void Ar_moleculardynamics::MD_initVel<UseAVX::True>()
	{
		auto const v = std::sqrt(3.0 * Tg_);

		myrandom::MyRand mr(-1.0, 1.0);

		for (auto n = 0; n < NumAtom_; n++) {
			double rndX = mr.myrand();
			double rndY = mr.myrand();
			double rndZ = mr.myrand();
			double tmp = 1.0 / std::sqrt(norm2(rndX, rndY, rndZ));
			rndX *= tmp;
			rndY *= tmp;
			rndZ *= tmp;

			// 方向はランダムに与える
			V_[n] = Ar_moleculardynamics::pack_t(v * rndX, v * rndY, v * rndZ, 0.0);
		}

		auto sx = 0.0;
		auto sy = 0.0;
		auto sz = 0.0;

		for (auto n = 0; n < NumAtom_; n++) {
			sx += V_[n][0];
			sy += V_[n][1];
			sz += V_[n][2];
		}

		sx /= static_cast<double>(NumAtom_);
		sy /= static_cast<double>(NumAtom_);
		sz /= static_cast<double>(NumAtom_);

		// 重心の並進運動を避けるために、速度の和がゼロになるように補正
		for (auto n = 0; n < NumAtom_; n++) {
			V_[n][0] -= sx;
			V_[n][1] -= sy;
			V_[n][2] -= sz;
		}
	}

	template<>
	inline void Ar_moleculardynamics::MD_initVel<UseAVX::False>()
	{
		auto const v = std::sqrt(3.0 * Tg_);

		myrandom::MyRand mr(-1.0, 1.0);

		for (auto n = 0; n < NumAtom_; n++) {
			double rndX = mr.myrand();
			double rndY = mr.myrand();
			double rndZ = mr.myrand();
			double tmp = 1.0 / std::sqrt(norm2(rndX, rndY, rndZ));
			rndX *= tmp;
			rndY *= tmp;
			rndZ *= tmp;

			// 方向はランダムに与える
			VX_[n] = v * rndX;
			VY_[n] = v * rndY;
			VZ_[n] = v * rndZ;
		}

		auto sx = 0.0;
		auto sy = 0.0;
		auto sz = 0.0;

		for (auto n = 0; n < NumAtom_; n++) {
			sx += VX_[n];
			sy += VY_[n];
			sz += VZ_[n];
		}

		sx /= static_cast<double>(NumAtom_);
		sy /= static_cast<double>(NumAtom_);
		sz /= static_cast<double>(NumAtom_);

		// 重心の並進運動を避けるために、速度の和がゼロになるように補正
		for (auto n = 0; n < NumAtom_; n++) {
			VX_[n] -= sx;
			VY_[n] -= sy;
			VZ_[n] -= sz;
		}
	}

	// #endregion privateメンバ関数 
}

#endif      // _AR_MOLECULARDYNAMICS_H_
