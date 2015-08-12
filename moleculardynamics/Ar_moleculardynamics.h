#ifdef _AR_MOLECULARDYNAMICS_H_
#define _AR_MOLECULARDYNAMICS_H_

#include <array>    // for std::array

namespace moleculardynamics {
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
        Ar_moleculardynamics() = default;

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~Ar_moleculardynamics() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数
        
        //! A public member function.
        /*!
            原子に働く力を計算する
        */
        void Calc_Forces();

        //! A public member function.
        /*!
            原子を移動させる
        */
        void Move_Atoms();
        //void Output_Data();

        // #endregion メンバ関数

        // #region メンバ変数

        //! A private member variable.
        /*!
            時間刻み
        */
        double dt;
        double const Tg;
        double const rc;
        double const rc2;
        double const rcm12;
        double const rcm6;
        double const Vrc;
        double const scale;

        double const dt2;
        double const alpha;

        double lat;
        double rho;
        double Tc;
        double Up;
        double Uk;
        double Utot;
        double t;

        const int ncp;
        const int Nc;

        int NumAtom;
        uint MD_iter;

        boost::array<double, ASIZE1> X;
        boost::array<double, ASIZE1> Y;
        boost::array<double, ASIZE1> Z;
        boost::array<double, ASIZE1> X1;
        boost::array<double, ASIZE1> Y1;
        boost::array<double, ASIZE1> Z1;
        boost::array<double, ASIZE1> FX;
        boost::array<double, ASIZE1> FY;
        boost::array<double, ASIZE1> FZ;
        boost::array<double, ASIZE1> VX;
        boost::array<double, ASIZE1> VY;
        boost::array<double, ASIZE1> VZ;

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        Ar_moleculardynamics() = delete;

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
