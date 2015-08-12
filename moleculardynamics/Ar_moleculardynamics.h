#ifdef _AR_MOLECULARDYNAMICS_H_
#define _AR_MOLECULARDYNAMICS_H_

#include <array>    // for std::array

namespace moleculardynamics {
    //! A class.
    /*!
        �A���S���ɑ΂��āA���q���͊w�V�~�����[�V�������s���N���X
    */
    class Ar_moleculardynamics final {
        // #region �R���X�g���N�^�E�f�X�g���N�^

    public:
        //! A constructor.
        /*!
            �R���X�g���N�^
        */
        Ar_moleculardynamics() = default;

        //! A destructor.
        /*!
            �f�t�H���g�f�X�g���N�^
        */
        ~Ar_moleculardynamics() = default;

        // #endregion �R���X�g���N�^�E�f�X�g���N�^

        // #region �����o�֐�
        
        //! A public member function.
        /*!
            ���q�ɓ����͂��v�Z����
        */
        void Calc_Forces();

        //! A public member function.
        /*!
            ���q���ړ�������
        */
        void Move_Atoms();
        //void Output_Data();

        // #endregion �����o�֐�

        // #region �����o�ϐ�

        //! A private member variable.
        /*!
            ���ԍ���
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

        // #region �֎~���ꂽ�R���X�g���N�^�E�����o�֐�

        //! A private constructor (deleted).
        /*!
            �f�t�H���g�R���X�g���N�^�i�֎~�j
        */
        Ar_moleculardynamics() = delete;

        //! A private copy constructor (deleted).
        /*!
            �R�s�[�R���X�g���N�^�i�֎~�j
        */
        Ar_moleculardynamics(Ar_moleculardynamics const &) = delete;

        //! A private member function (deleted).
        /*!
        operator=()�̐錾�i�֎~�j
        \param �R�s�[���̃I�u�W�F�N�g�i���g�p�j
        \return �R�s�[���̃I�u�W�F�N�g
        */
        Ar_moleculardynamics & operator=(Ar_moleculardynamics const &) = delete;

        // #endregion �֎~���ꂽ�R���X�g���N�^�E�����o�֐�
    };
}

#endif      // _AR_MOLECULARDYNAMICS_H_
