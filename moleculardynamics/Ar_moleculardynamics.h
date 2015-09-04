/*! \file Ar_moleculardynamics.h
	\brief アルゴンに対して、分子動力学シミュレーションを行うクラスの宣言

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#ifndef _AR_MOLECULARDYNAMICS_H_
#define _AR_MOLECULARDYNAMICS_H_

#pragma once

#include "../utility/property.h"
#include <array>					// for std::array
#include <cstdint>					// for std::int32_t
#include <vector>					// for std::vector

namespace moleculardynamics {
	using namespace utility;

    //! A class.
    /*!
        アルゴンに対して、分子動力学シミュレーションを行うクラス
    */
    class Ar_moleculardynamics final {
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

        //! A public member function.
        /*!
            原子に働く力を計算する
        */
        void Calc_Forces();

        //! A public member function.
        /*!
            n番目の原子に働く力を求める
        */
        float getForce(std::int32_t n) const;

        //! A public member function.
        /*!
            原子を移動させる
        */
        void Move_Atoms();

        //! A public member function.
        /*!
            再計算する
        */
        void recalc();

        // #endregion publicメンバ関数

		// #region privateメンバ関数

	private:
		//! A private member function.
		/*!
			原子の初期位置を決める
		*/
		void MD_initPos();

		//! A private member function.
		/*!
			原子の初期速度を決める
		*/
		void MD_initVel();

		//! A private member function.
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

	private:
		//! A private member variable (constant).
		/*!
			Woodcockの温度スケーリングの係数
		*/
		double const alpha = 0.2;

        //! A private member variable (constant).
        /*!
            時間刻み
        */
        double const dt = 0.001;

		//! A private member variable (constant).
		/*!
			時間刻みの二乗
		*/
		double const dt2;

		//! A private member variable.
		/*!
			格子定数
		*/
		double lat;

		//! A private member variable (constant).
		/*!
			スーパーセルの個数
		*/
		std::int32_t const Nc = 5;

		//! A private member variable.
		/*!
			n個目の原子に働く力のx成分
		*/
		std::vector<double> FX;

		//! A private member variable.
		/*!
			n個目の原子に働く力のy成分
		*/
		std::vector<double> FY;

		//! A private member variable.
		/*!
			n個目の原子に働く力のz成分
		*/
		std::vector<double> FZ;

		//! A private member variable.
		/*!
			MDのステップ数
		*/
		std::int32_t MD_iter = 1;

		//! A private member variable (constant).
		/*!
			相互作用を計算するセルの個数
		*/
		std::int32_t const ncp = 2;
		
		//! A private member variable.
		/*!
			原子の個数
		*/
		std::int32_t NumAtom;

		//! A private member variable (constant).
		/*!
			カットオフ半径
		*/
		double const rc = 2.5;

		//! A private member variable (constant).
		/*!
			カットオフ半径の2乗
		*/
		double const rc2;

		//! A private member variable (constant).
		/*!
			カットオフ半径の逆数の6乗
		*/
		double const rcm6;

		//! A private member variable (constant).
		/*!
			カットオフ半径の逆数の12乗
		*/
		double const rcm12;
		
		//! A private member variable.
		/*!
			密度
		*/
		double rho;

		//! A private member variable (constant).
		/*!
			格子定数のスケーリングの定数
		*/
		double const scale = 1.0;

		//! A private member variable.
		/*!
			時間	
		*/
		double t = 0.0;

		//! A private member variable.
		/*!
			計算された温度Tcalc
		*/
        double Tc;

		//! A private member variable.
		/*!
			与える温度Tgiven
		*/
        double Tg = 0.5;

		//! A private member variable.
		/*!
			運動エネルギー
		*/
		double Uk;
		
		//! A private member variable.
		/*!
			ポテンシャルエネルギー
		*/
		double Up;

		//! A private member variable.
		/*!
			総エネルギー
		*/
		double Utot;
		
		//! A private member variable (constant).
		/*!
			ポテンシャルエネルギーの打ち切り
		*/
        double const Vrc;
        
		//! A private member variable.
		/*!
			n個目の原子の速度のx成分
		*/
		std::vector<double> VX;

		//! A private member variable.
		/*!
			n個目の原子の速度のy成分
		*/
		std::vector<double> VY;

		//! A private member variable.
		/*!
			n個目の原子の速度のz成分
		*/
		std::vector<double> VZ;
		        
		//! A private member variable.
		/*!
			n個目の原子のx座標
		*/
        std::vector<double> X_;

		//! A private member variable.
		/*!
			n個目の原子の初期x座標
		*/
		std::vector<double> X1;
		
		//! A private member variable.
		/*!
			n個目の原子のy座標
		*/
		std::vector<double> Y_;

		//! A private member variable.
		/*!
			n個目の原子の初期y座標
		*/
		std::vector<double> Y1;
		
		//! A private member variable.
		/*!
			n個目の原子のz座標
		*/
		std::vector<double> Z_;
        
		//! A private member variable.
		/*!
			n個目の原子の初期z座標
		*/
		std::vector<double> Z1;

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
}

#endif      // _AR_MOLECULARDYNAMICS_H_
