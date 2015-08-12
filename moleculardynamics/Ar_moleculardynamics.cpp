#include "Ar_moleculardynamics.h"

namespace moleculardynamics {

MyMD_Sim::MyMDSimImpl::MyMDSimImpl()
 :	LogicalCore(GetLogicalProcessNumber()),
	Tg(2.0), rc(2.5), rc2(rc * rc),
	rcm12(std::pow(rc, -12.0)), rcm6(std::pow(rc, -6.0)),
	Vrc(4.0 * (rcm12 - rcm6)), scale(1.0),
	dt(0.001), dt2(dt * dt), alpha(0.2),
	ncp(2), Nc(5),
	ofsdata("M431_report_13.xyz"),
	ofstemp("M431_rep13_temp.txt")
{
	if (!ofsdata.is_open()) {
		throw std::domain_error("File open error!!\n");
	}
	if (!ofstemp.is_open()) {
		throw std::domain_error("File open error!!\n");
	}

    // initalize parameters
	lat = std::pow(2.0, 2.0 / 3.0) * scale;
    double uvol = lat * lat * lat; 
    rho = 4.0 / uvol;
    
	MD_iter = 1;

	t = 0.0;

	Uk = 0.0;
	Up = 0.0;
	Utot = 0.0;

	Tc = 0.0;

	init_MyMD_pos();
	init_MyMD_vel();

	lat *= boost::numeric_cast<double>(Nc);
}

MyMD_Sim::MyMDSimImpl::~MyMDSimImpl()
{
}

MyMD_Sim::MyMD_Sim() : pimpl_(new MyMDSimImpl)
{
}

MyMD_Sim::~MyMD_Sim()
{
}

void MyMD_Sim::MyMDSimImpl::init_MyMD_pos()
{
#ifdef _OPENMP
	omp_set_num_threads(LogicalCore);
#pragma omp parallel
#endif
	;

	double sx, sy, sz;
	int n = 0;

    // initial structure
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < Nc; i++) {
        for (int j = 0; j < Nc; j++) {
            for (int k = 0; k < Nc; k++) {
                // 基本セルをコピーする
                sx = boost::numeric_cast<double>(i) * lat; 
                sy = boost::numeric_cast<double>(j) * lat;
                sz = boost::numeric_cast<double>(k) * lat;

                // 基本セル内には4つの原子がある
                X[n] = sx;
                Y[n] = sy;
                Z[n] = sz;
                n++;

                X[n] = 0.5 * lat + sx;
                Y[n] = 0.5 * lat + sy;
                Z[n] = sz;
                n++;

                X[n] = sx;
                Y[n] = 0.5 * lat + sy;
                Z[n] = 0.5 * lat + sz;
                n++;

                X[n] = 0.5 * lat + sx; 
                Y[n] = sy;
                Z[n] = 0.5 * lat + sz;
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

    for (int n = 0; n < NumAtom; n++) {
        sx += X[n];
        sy += Y[n];
        sz += Z[n];
    }

    sx /= boost::numeric_cast<double>(NumAtom);
    sy /= boost::numeric_cast<double>(NumAtom);
    sz /= boost::numeric_cast<double>(NumAtom);

#ifdef _OPENMP
#pragma omp parallel for shared(sx, sy, sz)
#endif
    for (int n = 0; n < NumAtom; n++) {
        X[n] -= sx; 
        Y[n] -= sy; 
        Z[n] -= sz; 
    }
}

void MyMD_Sim::MyMDSimImpl::init_MyMD_vel()
{
#ifdef _OPENMP
	omp_set_num_threads(LogicalCore);
#pragma omp parallel
#endif
	;

	const boost::mt19937 randEngine(boost::numeric_cast<boost::uint32_t>(std::time(NULL)));
	boost::variate_generator<boost::mt19937, boost::uniform_real<> >
				myrand(randEngine, boost::uniform_real<>(-1.0, 1.0));

	const double v = std::sqrt(3.0 * Tg);

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int n = 0; n < NumAtom; n++) {
		double rndX = myrand();
		double rndY = myrand();
		double rndZ = myrand();
		double tmp = 1.0 / std::sqrt(MyMD_Sim::norm2(rndX, rndY, rndZ));
		rndX *= tmp;
		rndY *= tmp;
		rndZ *= tmp;

		// 方向はランダムに与える
		VX[n] = v * rndX;
		VY[n] = v * rndY;
		VZ[n] = v * rndZ;
	}

	double sx = 0.0;
	double sy = 0.0;
	double sz = 0.0;

	for (int n = 0; n < NumAtom; n++) {
		sx += VX[n];
		sy += VY[n];
		sz += VZ[n];
	}

	sx /= boost::numeric_cast<double>(NumAtom);
	sy /= boost::numeric_cast<double>(NumAtom);
	sz /= boost::numeric_cast<double>(NumAtom);

	// 重心の並進運動を避けるために、速度の和がゼロになるように補正
#ifdef _OPENMP
#pragma omp parallel for shared(sx, sy, sz)
#endif
	for (int n = 0; n < NumAtom; n++) {
		VX[n] -= sx;
		VY[n] -= sy;
		VZ[n] -= sz;
	}
}

double MyMD_Sim::norm2(const double x, const double y, const double z)
{
	return (x * x + y * y + z * z);
}

void MyMD_Sim::Calc_Forces()
{
	for (uint n = 0; n < pimpl_->NumAtom; n++) {
		for (uint m = 0; m < pimpl_->NumAtom; m++) {

			// ±ncp分のセル内の原子との相互作用を計算
			for (int i = - pimpl_->ncp; i <= pimpl_->ncp; i++) {
				for (int j = - pimpl_->ncp; j <= pimpl_->ncp; j++) {
					for (int k = - pimpl_->ncp; k <= pimpl_->ncp; k++) {
						double sx = boost::numeric_cast<double>(i) * pimpl_->lat; 
						double sy = boost::numeric_cast<double>(j) * pimpl_->lat;
						double sz = boost::numeric_cast<double>(k) * pimpl_->lat;

						// 自分自身との相互作用を排除
						if (n != m || i != 0 || j != 0 || k != 0) {
							double dx = pimpl_->X[n] - (pimpl_->X[m] + sx);
							double dy = pimpl_->Y[n] - (pimpl_->Y[m] + sy);
							double dz = pimpl_->Z[n] - (pimpl_->Z[m] + sz);

							double r2 = MyMD_Sim::norm2(dx, dy, dz);
							// 打ち切り距離内であれば計算
							if (r2 <= pimpl_->rc2) {
								double r = sqrt(r2);
								double rm6 = 1.0 / (r2 * r2 * r2);
								double rm7 = rm6 / r;
								double rm12 = rm6 * rm6;
								double rm13 = rm12 / r;

								double Fr = 48.0 * rm13 - 24.0 * rm7;

								pimpl_->FX[n] += dx / r * Fr;
								pimpl_->FY[n] += dy / r * Fr;
								pimpl_->FZ[n] += dz / r * Fr;

								// 0.5 for double counting
								// エネルギーの計算、ただし二重計算のために0.5をかけておく
								pimpl_->Up += 0.5 * (4.0 * (rm12 - rm6) - pimpl_->Vrc);
							}							
						}
					}
				}
			}
		}
	}
#endif
}

void MyMD_Sim::Move_Atoms()
{
#ifdef _OPENMP
	omp_set_num_threads(pimpl_->LogicalCore);
#pragma omp parallel
#endif
	;

	// calculate temperture
	pimpl_->Uk = 0.0;

	for (int n = 0; n < pimpl_->NumAtom; n++) {
		double v2 = MyMD_Sim::norm2(pimpl_->VX[n], pimpl_->VY[n], pimpl_->VZ[n]);
		pimpl_->Uk += v2;
	}

	// 運動エネルギーの計算
	pimpl_->Uk *= 0.5;
	// 全エネルギーの計算
	pimpl_->Utot = pimpl_->Up + pimpl_->Uk;
	// 温度の計算
	pimpl_->Tc = pimpl_->Uk / (1.5 * boost::numeric_cast<double>(pimpl_->NumAtom));

	double s = std::sqrt((pimpl_->Tg + pimpl_->alpha * (pimpl_->Tc - pimpl_->Tg)) / pimpl_->Tc);

	if (pimpl_->MD_iter == 1) {
		// update the coordinates by the second order Euler method
		// 最初のステップだけ修正Euler法で時間発展
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int n = 0; n < pimpl_->NumAtom; n++) {
			pimpl_->X1[n] = pimpl_->X[n];
			pimpl_->Y1[n] = pimpl_->Y[n];
			pimpl_->Z1[n] = pimpl_->Z[n];

			// scaling of velocity
			pimpl_->VX[n] *= s;
			pimpl_->VY[n] *= s;
			pimpl_->VZ[n] *= s;

			// update coordinates and velocity
			pimpl_->X[n] += pimpl_->dt * pimpl_->VX[n] + 0.5 * pimpl_->FX[n] * pimpl_->dt2;
			pimpl_->Y[n] += pimpl_->dt * pimpl_->VY[n] + 0.5 * pimpl_->FY[n] * pimpl_->dt2;
			pimpl_->Z[n] += pimpl_->dt * pimpl_->VZ[n] + 0.5 * pimpl_->FZ[n] * pimpl_->dt2;

			pimpl_->VX[n] += pimpl_->dt * pimpl_->FX[n];
			pimpl_->VY[n] += pimpl_->dt * pimpl_->FY[n];
			pimpl_->VZ[n] += pimpl_->dt * pimpl_->FZ[n];
		}
	} else {
		// update the coordinates by the Verlet method
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int n = 0; n < pimpl_->NumAtom; n++) {
			double rtmp[3];
			rtmp[0] = pimpl_->X[n];
			rtmp[1] = pimpl_->Y[n];
			rtmp[2] = pimpl_->Z[n];

			// update coordinates and velocity
			// Verlet法の座標更新式において速度成分を抜き出し、その部分をスケールする
#ifdef NVT
			pimpl_->X[n] += s * (pimpl_->X[n] - pimpl_->X1[n]) + pimpl_->FX[n] * pimpl_->dt2;
			pimpl_->Y[n] += s * (pimpl_->Y[n] - pimpl_->Y1[n]) + pimpl_->FY[n] * pimpl_->dt2;
			pimpl_->Z[n] += s * (pimpl_->Z[n] - pimpl_->Z1[n]) + pimpl_->FZ[n] * pimpl_->dt2;
#elif NVE
			pimpl_->X[n] = 2.0 * pimpl_->X[n] - pimpl_->X1[n] + pimpl_->FX[n] * pimpl_->dt2;
			pimpl_->Y[n] = 2.0 * pimpl_->Y[n] - pimpl_->Y1[n] + pimpl_->FY[n] * pimpl_->dt2;
			pimpl_->Z[n] = 2.0 * pimpl_->Z[n] - pimpl_->Z1[n] + pimpl_->FZ[n] * pimpl_->dt2;
#endif
			pimpl_->VX[n] = 0.5 * (pimpl_->X[n] - pimpl_->X1[n]) / pimpl_->dt;
			pimpl_->VY[n] = 0.5 * (pimpl_->Y[n] - pimpl_->Y1[n]) / pimpl_->dt;
			pimpl_->VZ[n] = 0.5 * (pimpl_->Z[n] - pimpl_->Z1[n]) / pimpl_->dt;

			pimpl_->X1[n] = rtmp[0];
			pimpl_->Y1[n] = rtmp[1];
			pimpl_->Z1[n] = rtmp[2];
		}
	}

	// consider the periodic boundary condination
	// セルの外側に出たら座標をセル内に戻す
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int n = 0; n < pimpl_->NumAtom; n++) {
		if (pimpl_->X[n] > pimpl_->lat) {
			pimpl_->X[n] -= pimpl_->lat;
			pimpl_->X1[n] -= pimpl_->lat;
		} else if (pimpl_->X[n] < 0.0) {
			pimpl_->X[n] += pimpl_->lat;
			pimpl_->X1[n] += pimpl_->lat;
		}
		if (pimpl_->Y[n] > pimpl_->lat) {
			pimpl_->Y[n] -= pimpl_->lat;
			pimpl_->Y1[n] -= pimpl_->lat;
		} else if (pimpl_->Y[n] < 0.0) {
			pimpl_->Y[n] += pimpl_->lat;
			pimpl_->Y1[n] += pimpl_->lat;
		}
		if (pimpl_->Z[n] > pimpl_->lat) {
			pimpl_->Z[n] -= pimpl_->lat;
			pimpl_->Z1[n] -= pimpl_->lat;
		} else if (pimpl_->Z[n] < 0.0) {
			pimpl_->Z[n] += pimpl_->lat;
			pimpl_->Z1[n] += pimpl_->lat;
		}
	}

	// 繰り返し回数と時間を増加
	pimpl_->t = boost::numeric_cast<double>(pimpl_->MD_iter) * pimpl_->dt;
	pimpl_->MD_iter++;
}

void MyMD_Sim::Output_Data()
{
	pimpl_->ofsdata << pimpl_->NumAtom << std::endl;
	pimpl_->ofsdata << boost::format("iter=%02d ") % pimpl_->MD_iter
				   << boost::format("t=%8.3f ") % pimpl_->t
				   << boost::format("Uk=%8.5f ") % pimpl_->Uk
				   << boost::format("Up=%8.5f ") % pimpl_->Up
				   << boost::format("Utot=%8.5f ") % pimpl_->Utot
				   << boost::format("Tg=%6.3f ") % pimpl_->Tg
				   << boost::format("Tc=%6.3f") % pimpl_->Tc
				   << std::endl;
	
	for (int n = 0; n < pimpl_->NumAtom; n++) {
		pimpl_->ofsdata << "Ar   "
					   << boost::format("%8.5f ") % (pimpl_->X[n] * 3.0)
					   << boost::format("%8.5f ") % (pimpl_->Y[n] * 3.0)
					   << boost::format("%8.5f") % (pimpl_->Z[n] * 3.0)
					   << "   "
					   << boost::format("%8.5f ") % pimpl_->VX[n]
					   << boost::format("%8.5f ") % pimpl_->VY[n]
					   << boost::format("%8.5f") % pimpl_->VZ[n]
					   << std::endl;
	}

	pimpl_->ofstemp << boost::format("%d\t") % pimpl_->MD_iter
					  << boost::format("%8.5f\t") % pimpl_->Tg
					  << boost::format("%8.5f") % pimpl_->Tc
					  << std::endl;
}
