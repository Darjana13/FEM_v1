#include "Gauss.h"

double mult(vector<double>& a, vector<double>& b)
{
	int n = a.size();
	double res = 0;
	for (int i = 0; i < n; i++)
	{
		res += a[i] * b[i];
	}
	return res;
}

double Integrate_Gauss3Method2D::IntegrateFunc(Node2D& from, Node2D& to, func_2coord_to_1& func)
{
	Node2D ksi_from, ksi_to;

	double ksi_x, ksi_y;
	double x_centre = (to.x + from.x) / 2, y_centre = (to.y + from.y) / 2;
	double h_x = (to.x - from.x) / 2, h_y = (to.y - from.y) / 2;

	int p_id = 0;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		{
			intagrate_points_2D[p_id].coord[0] = x_centre + h_x * basic_point[i];
			intagrate_points_2D[p_id].coord[1] = y_centre + h_y * basic_point[j];
			intagrate_points_2D[p_id].koef = basic_coef[i] * basic_coef[j];
			p_id++;
		}

	double res = 0;
	for (int p_id = 0; p_id < 3 * 3; p_id++)
	{
		res += h_x * h_y * intagrate_points_2D[p_id].koef * func(intagrate_points_2D[p_id].coord[0], intagrate_points_2D[p_id].coord[1]);
	}
	return res;
}

void Integrate_Gauss3Method2D::Init()
{
	phi_in_points.resize(2);
	derivative_phi_in_points.resize(2);
	double ksi_x, ksi_y;
	intagrate_points_2D.resize(3 * 3);
	int p_id = 0;
	for (int i = 0; i < 3; i++)
			for (int k = 0; k < 3; k++)
			{
				intagrate_points_2D[p_id].coord[0] = basic_point[i];
				intagrate_points_2D[p_id].coord[1] = basic_point[k];
				intagrate_points_2D[p_id].koef = basic_coef[i] * basic_coef[k];
				p_id++;
			}
}


//void Integrate_Gauss3Method3D::Init(vector<basis_func>& basis, vector<basis_func>& derivative_basis)
//{
//	double ksi_x, ksi_y, ksi_z;
//
//	int p_id = 0;
//	for (int i = 0; i < 3; i++)
//		for (int j = 0; j < 3; j++)
//			for (int k = 0; k < 3; k++)
//			{
//				intagrate_points_3D[p_id].coord[0] = basic_point[i];
//				intagrate_points_3D[p_id].coord[1] = basic_point[j];
//				intagrate_points_3D[p_id].coord[2] = basic_point[k];
//				intagrate_points_3D[p_id].koef = basic_coef[i] * basic_coef[j] * basic_coef[k];
//				ksi_x = basic_point[i];
//				ksi_y = basic_point[j];
//				ksi_z = basic_point[k];
//				for (int node_id = 0; node_id < 8; node_id++)
//				{
//					intagrate_points_3D[p_id].derivative[node_id][0] = derivative_basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//					intagrate_points_3D[p_id].derivative[node_id][1] = basis[mu(node_id)](ksi_x) * derivative_basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//					intagrate_points_3D[p_id].derivative[node_id][2] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * derivative_basis[teta(node_id)](ksi_z);
//					intagrate_points_3D[p_id].basis_func[node_id] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//				}
//				p_id++;
//			}
//
//	double ksi, eta;
//	
//	for (int face_id = 0; face_id < 6; face_id++)
//	{
//		p_id = 0;
//		switch (face_id)
//		{
//		case 0:
//			for (int i = 0; i < 3; i++)
//				for (int j = 0; j < 3; j++)
//				{
//					intagrate_points_face[face_id][p_id].coord[0] = -1;
//					intagrate_points_face[face_id][p_id].coord[1] = basic_point[i];
//					intagrate_points_face[face_id][p_id].coord[2] = basic_point[j];
//					intagrate_points_face[face_id][p_id].koef = basic_coef[i] * basic_coef[j]; 
//					ksi_x = -1;
//					ksi_y = basic_point[i];
//					ksi_z = basic_point[j];
//					for (int node_id = 0; node_id < 8; node_id++) // ??? зачем 8
//					{
//						intagrate_points_face[face_id][p_id].derivative[node_id][0] = derivative_basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].derivative[node_id][1] = basis[mu(node_id)](ksi_x) * derivative_basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].derivative[node_id][2] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * derivative_basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].basis_func[node_id] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//					}
//					p_id++;
//				}
//			break;
//		case 1:
//			for (int i = 0; i < 3; i++)
//				for (int j = 0; j < 3; j++)
//				{
//					intagrate_points_face[face_id][p_id].coord[0] = 1;
//					intagrate_points_face[face_id][p_id].coord[1] = basic_point[i];
//					intagrate_points_face[face_id][p_id].coord[2] = basic_point[j];
//					intagrate_points_face[face_id][p_id].koef = basic_coef[i] * basic_coef[j];
//					ksi_x = 1;
//					ksi_y = basic_point[i];
//					ksi_z = basic_point[j];
//					for (int node_id = 0; node_id < 8; node_id++)
//					{
//						intagrate_points_face[face_id][p_id].derivative[node_id][0] = derivative_basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].derivative[node_id][1] = basis[mu(node_id)](ksi_x) * derivative_basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].derivative[node_id][2] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * derivative_basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].basis_func[node_id] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//					}
//					p_id++;
//				}
//			break;
//		case 2:
//			for (int i = 0; i < 3; i++)
//				for (int j = 0; j < 3; j++)
//				{
//					intagrate_points_face[face_id][p_id].coord[0] = basic_point[i];
//					intagrate_points_face[face_id][p_id].coord[1] = -1;
//					intagrate_points_face[face_id][p_id].coord[2] = basic_point[j];
//					intagrate_points_face[face_id][p_id].koef = basic_coef[i] * basic_coef[j];
//					ksi_x = basic_point[i];
//					ksi_y = -1;
//					ksi_z = basic_point[j];
//					for (int node_id = 0; node_id < 8; node_id++)
//					{
//						intagrate_points_face[face_id][p_id].derivative[node_id][0] = derivative_basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].derivative[node_id][1] = basis[mu(node_id)](ksi_x) * derivative_basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].derivative[node_id][2] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * derivative_basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].basis_func[node_id] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//					}
//					p_id++;
//				}
//			break;
//		case 3:
//			for (int i = 0; i < 3; i++)
//				for (int j = 0; j < 3; j++)
//				{
//					intagrate_points_face[face_id][p_id].coord[0] = basic_point[i];
//					intagrate_points_face[face_id][p_id].coord[1] = 1;
//					intagrate_points_face[face_id][p_id].coord[2] = basic_point[j];
//					intagrate_points_face[face_id][p_id].koef = basic_coef[i] * basic_coef[j];
//					ksi_x = basic_point[i];
//					ksi_y = 1;
//					ksi_z = basic_point[j];
//					for (int node_id = 0; node_id < 8; node_id++)
//					{
//						intagrate_points_face[face_id][p_id].derivative[node_id][0] = derivative_basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].derivative[node_id][1] = basis[mu(node_id)](ksi_x) * derivative_basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].derivative[node_id][2] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * derivative_basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].basis_func[node_id] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//					}
//					p_id++;
//				}
//			break;
//		case 4:
//			for (int i = 0; i < 3; i++)
//				for (int j = 0; j < 3; j++)
//				{
//					intagrate_points_face[face_id][p_id].coord[0] = basic_point[i];
//					intagrate_points_face[face_id][p_id].coord[1] = basic_point[j];
//					intagrate_points_face[face_id][p_id].coord[2] = -1;
//					intagrate_points_face[face_id][p_id].koef = basic_coef[i] * basic_coef[j];
//					ksi_x = basic_point[i];
//					ksi_y = basic_point[j];
//					ksi_z = -1;
//					for (int node_id = 0; node_id < 8; node_id++)
//					{
//						intagrate_points_face[face_id][p_id].derivative[node_id][0] = derivative_basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].derivative[node_id][1] = basis[mu(node_id)](ksi_x) * derivative_basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].derivative[node_id][2] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * derivative_basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].basis_func[node_id] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//					}
//					p_id++;
//				}
//			break;
//		case 5:
//			for (int i = 0; i < 3; i++)
//				for (int j = 0; j < 3; j++)
//				{
//					intagrate_points_face[face_id][p_id].coord[0] = basic_point[i];
//					intagrate_points_face[face_id][p_id].coord[1] = basic_point[j];
//					intagrate_points_face[face_id][p_id].coord[2] = 1;
//					intagrate_points_face[face_id][p_id].koef = basic_coef[i] * basic_coef[j];
//					ksi_x = basic_point[i];
//					ksi_y = basic_point[j];
//					ksi_z = 1;
//					for (int node_id = 0; node_id < 8; node_id++)
//					{
//						intagrate_points_face[face_id][p_id].derivative[node_id][0] = derivative_basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].derivative[node_id][1] = basis[mu(node_id)](ksi_x) * derivative_basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].derivative[node_id][2] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * derivative_basis[teta(node_id)](ksi_z);
//						intagrate_points_face[face_id][p_id].basis_func[node_id] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//					}
//					p_id++;
//				}
//			break;
//		default:
//			break;
//		}
//	}
//}
//
//double Integrate_Gauss3Method3D::CalcdetJ(int p_id, vector<Node3D> &element)
//{
//	vector<vector<double>> J(3, vector<double>(3));
//	for (int i = 0; i < 8; i++)
//	{
//		J[0][0] += element[i].x * intagrate_points_3D[p_id].derivative[i][0]; //  d_xi[i];
//		J[0][1] += element[i].x * intagrate_points_3D[p_id].derivative[i][1]; //  d_eta[i];
//		J[0][2] += element[i].x * intagrate_points_3D[p_id].derivative[i][2]; //  d_zeta[i];
//
//		J[1][0] += element[i].y * intagrate_points_3D[p_id].derivative[i][0]; //  d_xi[i];
//		J[1][1] += element[i].y * intagrate_points_3D[p_id].derivative[i][1]; //  d_eta[i];
//		J[1][2] += element[i].y * intagrate_points_3D[p_id].derivative[i][2]; //  d_zeta[i];
//
//		J[2][0] += element[i].z * intagrate_points_3D[p_id].derivative[i][0]; //  d_xi[i];
//		J[2][1] += element[i].z * intagrate_points_3D[p_id].derivative[i][1]; //  d_eta[i];
//		J[2][2] += element[i].z * intagrate_points_3D[p_id].derivative[i][2]; //  d_zeta[i];
//	}
//
//	// вычисл€ем якобиан (определитель)
//	double det_J = J[0][0] * J[1][1] * J[2][2] - J[0][0] * J[1][2] * J[2][1] + J[1][0] * J[2][1] * J[0][2]
//		- J[1][0] * J[0][1] * J[2][2] + J[2][0] * J[0][1] * J[1][2] - J[2][0] * J[1][1] * J[0][2];
//
//	return det_J;
//}
//
//double Integrate_Gauss3Method3D::CalcdetJ_face(int p_id, vector<Node3D>& face, int face_id)
//{
//	vector<double> dcoord_dksi(3), dcoord_dnu(3);
//	double J_11 = 0, J_12 = 0, J_22 = 0;
//	for (int node_i = 0; node_i < 4; node_i++)
//	{
//		int node_3D = face_to_node[face_id][node_i];
//		dcoord_dksi[0] += face[node_i].x * intagrate_points_face[face_id][p_id].derivative[node_3D][face_to_derivative[face_id][0]]; //  d_xi[i];
//		dcoord_dnu[0]  += face[node_i].x * intagrate_points_face[face_id][p_id].derivative[node_3D][face_to_derivative[face_id][1]];  //  d_eta[i];
//		dcoord_dksi[1] += face[node_i].y * intagrate_points_face[face_id][p_id].derivative[node_3D][face_to_derivative[face_id][0]]; //  d_xi[i];
//		dcoord_dnu[1]  += face[node_i].y * intagrate_points_face[face_id][p_id].derivative[node_3D][face_to_derivative[face_id][1]];  //  d_eta[i];
//		dcoord_dksi[2] += face[node_i].z * intagrate_points_face[face_id][p_id].derivative[node_3D][face_to_derivative[face_id][0]]; //  d_xi[i];
//		dcoord_dnu[2]  += face[node_i].z * intagrate_points_face[face_id][p_id].derivative[node_3D][face_to_derivative[face_id][1]];  //  d_eta[i];
//	}
//	for (int coord_i = 0; coord_i < 3; coord_i++)
//	{
//		J_11 += dcoord_dksi[coord_i] * dcoord_dksi[coord_i];
//		J_12 += dcoord_dksi[coord_i] * dcoord_dnu[coord_i];
//		J_22 += dcoord_dnu[coord_i] * dcoord_dnu[coord_i];
//	}
//
//	// вычисл€ем якобиан (определитель)
//	double det_J = J_11 * J_22 - J_12 * J_12;
//
//	return det_J;
//}
//
//double Integrate_Gauss3Method3D::CalcJ_1(vector<vector<double>>& J_1_T, int p_id, vector<Node3D> &element)
//{
//	vector<vector<double>> J(3, vector<double>(3));
//	for (int i = 0; i < 8; i++)
//	{
//		J[0][0] += element[i].x * intagrate_points_3D[p_id].derivative[i][0]; //  d_xi[i];
//		J[0][1] += element[i].x * intagrate_points_3D[p_id].derivative[i][1]; //  d_eta[i];
//		J[0][2] += element[i].x * intagrate_points_3D[p_id].derivative[i][2]; //  d_zeta[i];
//
//		J[1][0] += element[i].y * intagrate_points_3D[p_id].derivative[i][0]; //  d_xi[i];
//		J[1][1] += element[i].y * intagrate_points_3D[p_id].derivative[i][1]; //  d_eta[i];
//		J[1][2] += element[i].y * intagrate_points_3D[p_id].derivative[i][2]; //  d_zeta[i];
//
//		J[2][0] += element[i].z * intagrate_points_3D[p_id].derivative[i][0]; //  d_xi[i];
//		J[2][1] += element[i].z * intagrate_points_3D[p_id].derivative[i][1]; //  d_eta[i];
//		J[2][2] += element[i].z * intagrate_points_3D[p_id].derivative[i][2]; //  d_zeta[i];
//	}
//
//	// вычисл€ем якобиан (определитель)
//	double det_J = J[0][0] * J[1][1] * J[2][2] - J[0][0] * J[1][2] * J[2][1] + J[1][0] * J[2][1] * J[0][2]
//		- J[1][0] * J[0][1] * J[2][2] + J[2][0] * J[0][1] * J[1][2] - J[2][0] * J[1][1] * J[0][2];
//
//	// матрица, обратна€ к транспонированной матрице якоби
//	J_1_T[0][0] = (J[1][1] * J[2][2] - J[2][1] * J[1][2]) / det_J;
//	J_1_T[1][0] = (J[2][1] * J[0][2] - J[0][1] * J[2][2]) / det_J;
//	J_1_T[2][0] = (J[0][1] * J[1][2] - J[1][1] * J[0][2]) / det_J;
//	J_1_T[0][1] = (-J[1][0] * J[2][2] + J[2][0] * J[1][2]) / det_J;
//	J_1_T[1][1] = (J[0][0] * J[2][2] - J[2][0] * J[0][2]) / det_J;
//	J_1_T[2][1] = (-J[0][0] * J[1][2] + J[1][0] * J[0][2]) / det_J;
//	J_1_T[0][2] = (J[1][0] * J[2][1] - J[2][0] * J[1][1]) / det_J;
//	J_1_T[1][2] = (-J[0][0] * J[2][1] + J[2][0] * J[0][1]) / det_J;
//	J_1_T[2][2] = (J[0][0] * J[1][1] - J[1][0] * J[0][1]) / det_J;
//
//	return det_J;
//}
//
//int Integrate_Gauss3Method3D::GetLocG(vector<Node3D> &element, vector<vector<double>> &locG)
//{
//	for (int i = 0; i < 8; i++)
//	{
//		for (int j = 0; j < 8; j++)
//		{
//			locG[i][j] = 0;
//		}
//	}
//	for (int p_id = 0; p_id < 27; p_id++)
//	{
//		vector<vector<double>> J_1(3, vector<double>(3));
//		double det_J = CalcJ_1(J_1, p_id, element);
//
//		// intagrate_points_3D[p_id].derivative[i] - градиент i-ой базисной функции в точке p_id
//		vector<vector<double>> modifided_grad(8, vector<double>(3)); 
//		double modifided_coef = intagrate_points_3D[p_id].koef * abs(det_J);
//		for (int i = 0; i < 8; i++)
//		{
//			// modifided_grad[i] = J_1 * intagrate_points_3D[p_id].derivative[i]
//			for (int j = 0; j < 3; j++)
//			{
//				for (int k = 0; k < 3; k++)
//				{
//					modifided_grad[i][j] += J_1[j][k] * intagrate_points_3D[p_id].derivative[i][k];
//				}
//			}
//		}
//
//		for (int i = 0; i < 8; i++)
//		{
//			for (int j = 0; j < 8; j++)
//			{
//				locG[i][j] += mult(modifided_grad[i], modifided_grad[j]) * modifided_coef;
//			}
//		}
//
//	}
//
//	return 0;
//}
//
//int Integrate_Gauss3Method3D::GetLocM(vector<Node3D>& element, vector<vector<double>>& locM)
//{
//	for (int i = 0; i < 8; i++)
//	{
//		for (int j = 0; j < 8; j++)
//		{
//			locM[i][j] = 0;
//		}
//	}
//	for (int p_id = 0; p_id < 27; p_id++)
//	{
//		double det_J = CalcdetJ(p_id, element);
//
//		// intagrate_points_3D[p_id].basis_func[i] - i-а€ базисна€ функции в точке p_id
//		double modifided_coef = intagrate_points_3D[p_id].koef * abs(det_J);
//		for (int i = 0; i < 8; i++)
//		{
//			for (int j = 0; j < 8; j++)
//			{
//				locM[i][j] += intagrate_points_3D[p_id].basis_func[i] * intagrate_points_3D[p_id].basis_func[j] * modifided_coef;
//			}
//		}
//	}
//
//	return 0;
//}
//
//int Integrate_Gauss3Method3D::GetLocM_face(vector<Node3D>& face, vector<vector<double>>& locM, int face_id)
//{
//	if (locM.size() != 4)
//	{
//		cout << "error. Integrate_Gauss3Method3D::GetLocM_face vector<vector<double>>& locM must have size [4][4]" << endl;
//	}
//	for (int i = 0; i < 4; i++)
//	{
//		if (locM[i].size() != 4)
//		{
//			cout << "error. Integrate_Gauss3Method3D::GetLocM_face vector<vector<double>>& locM must have size [4][4]" << endl;
//		}
//		for (int j = 0; j < 4; j++)
//		{
//			locM[i][j] = 0;
//		}
//	}
//	for (int p_id = 0; p_id < 3*3; p_id++)
//	{
//		double det_J = CalcdetJ_face(p_id, face, face_id);
//		det_J = sqrt(abs(det_J)); // ??? почему корень
//		//cout << " det J " << det_J << endl;
//		// intagrate_points_3D[p_id].basis_func[i] - i-а€ базисна€ функции в точке p_id
//		double modifided_coef = intagrate_points_face[face_id][p_id].koef * abs(det_J);
//		
//		int i_3D, j_3D;
//		for (int i = 0; i < 4; i++)
//		{
//			i_3D = face_to_node[face_id][i];
//			for (int j = 0; j < 4; j++)
//			{
//				j_3D = face_to_node[face_id][j];
//				locM[i][j] += intagrate_points_face[face_id][p_id].basis_func[i_3D] * intagrate_points_face[face_id][p_id].basis_func[j_3D] * modifided_coef;
//			}
//		}
//	}
//
//	return 0;
//}
//
//
//double Integrate_Gauss3Method3D::IntegrateFunc(Node3D& from, Node3D& to, func_3coord_to_1& func)
//{
//	Node3D ksi_from, ksi_to;
//
//	double ksi_x, ksi_y, ksi_z;
//	double x_centre = (to.x + from.x) / 2, y_centre = (to.y + from.y) / 2, z_centre = (to.z + from.z) / 2;
//	double h_x = (to.x - from.x) / 2, h_y = (to.y - from.y) / 2, h_z = (to.z - from.z) / 2;
//
//	int p_id = 0;
//	for (int i = 0; i < 3; i++)
//		for (int j = 0; j < 3; j++)
//			for (int k = 0; k < 3; k++)
//			{
//				intagrate_points_3D[p_id].coord[0] = x_centre + h_x * basic_point[i];
//				intagrate_points_3D[p_id].coord[1] = y_centre + h_y * basic_point[j];
//				intagrate_points_3D[p_id].coord[2] = z_centre + h_z * basic_point[k];
//				intagrate_points_3D[p_id].koef = basic_coef[i] * basic_coef[j] * basic_coef[k];
//				/*ksi_x = getKsi(basic_point[i], -1, 1);
//				ksi_y = getKsi(basic_point[j], -1, 1);
//				ksi_z = getKsi(basic_point[k], -1, 1);
//				for (int node_id = 0; node_id < 8; node_id++)
//				{
//					intagrate_points_3D[p_id].derivative[node_id][0] = derivative_basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//					intagrate_points_3D[p_id].derivative[node_id][1] = basis[mu(node_id)](ksi_x) * derivative_basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//					intagrate_points_3D[p_id].derivative[node_id][2] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * derivative_basis[teta(node_id)](ksi_z);
//					intagrate_points_3D[p_id].basis_func[node_id] = basis[mu(node_id)](ksi_x) * basis[nu(node_id)](ksi_y) * basis[teta(node_id)](ksi_z);
//				}*/
//				p_id++;
//			}
//
//	double res = 0;
//	for (int p_id = 0; p_id < 3 * 3 * 3; p_id++)
//	{
//		res += h_x * h_y * h_z * intagrate_points_3D[p_id].koef * func(intagrate_points_3D[p_id].coord[0], intagrate_points_3D[p_id].coord[1], intagrate_points_3D[p_id].coord[2]);
//	}
//	return res;
//}
//
////----------------------------------------------------------------------------------------------------------
//
//void Integrate_Gauss3Method2D::Init(vector<basis_func>& basis, vector<basis_func>& derivative_basis)
//{
//	phi_in_points.resize(2);
//	derivative_phi_in_points.resize(2);
//	double ksi_x, ksi_y;
//	intagrate_points_2D.resize(3 * 3);
//	int p_id = 0;
//	for (int i = 0; i < 3; i++)
//			for (int k = 0; k < 3; k++)
//			{
//				intagrate_points_2D[p_id].coord[0] = basic_point[i];
//				intagrate_points_2D[p_id].coord[1] = basic_point[k];
//				intagrate_points_2D[p_id].koef = basic_coef[i] * basic_coef[k];
//				ksi_x = getKsi(basic_point[i], -1, 1);
//				ksi_y = getKsi(basic_point[k], -1, 1);
//				for (int node_id = 0; node_id < 4; node_id++)
//				{
//					intagrate_points_2D[p_id].derivative[node_id][0] = derivative_basis[node_id%2](ksi_x) * basis[node_id/2](ksi_y);
//					intagrate_points_2D[p_id].derivative[node_id][1] = basis[node_id % 2](ksi_x) * derivative_basis[node_id / 2](ksi_y);
//					intagrate_points_2D[p_id].basis_func[node_id] = basis[node_id % 2](ksi_x) * basis[node_id / 2](ksi_y);
//				}
//				p_id++;
//			}
//}
//
//double Integrate_Gauss3Method2D::CalcdetJ(int p_id, vector<Node2D>& element)
//{
//	vector<vector<double>> J(2, vector<double>(2));
//	for (int i = 0; i < 4; i++)
//	{
//		J[0][0] += element[i].x * intagrate_points_2D[p_id].derivative[i][0]; //  d_xi[i];
//		J[0][1] += element[i].x * intagrate_points_2D[p_id].derivative[i][1]; //  d_eta[i];
//		
//		J[1][0] += element[i].y * intagrate_points_2D[p_id].derivative[i][0]; //  d_xi[i];
//		J[1][1] += element[i].y * intagrate_points_2D[p_id].derivative[i][1]; //  d_eta[i];
//	}
//
//	// вычисл€ем якобиан (определитель)
//	double det_J = J[0][0] * J[1][1] - J[1][0] * J[0][1];
//
//	return det_J;
//}
//
//double Integrate_Gauss3Method2D::CalcJ_1(vector<vector<double>>& J_1_T, int p_id, vector<Node2D>& element)
//{
//	vector<vector<double>> J(2, vector<double>(2));
//	for (int i = 0; i < 8; i++)
//	{
//		J[0][0] += element[i].x * intagrate_points_2D[p_id].derivative[i][0]; //  d_xi[i];
//		J[0][1] += element[i].x * intagrate_points_2D[p_id].derivative[i][1]; //  d_eta[i];
//		
//		J[1][0] += element[i].y * intagrate_points_2D[p_id].derivative[i][0]; //  d_xi[i];
//		J[1][1] += element[i].y * intagrate_points_2D[p_id].derivative[i][1]; //  d_eta[i];
//	}
//
//	// вычисл€ем якобиан (определитель)
//	double det_J = J[0][0] * J[1][1] - J[1][0] * J[0][1];
//
//	// матрица, обратна€ к транспонированной матрице якоби
//	J_1_T[0][0] = J[1][1] / det_J;
//	J_1_T[1][0] = - J[0][1] / det_J;
//
//	J_1_T[0][1] = -J[1][0] / det_J;
//	J_1_T[1][1] = J[0][0]  / det_J;
//
//	return det_J;
//}
//
//int Integrate_Gauss3Method2D::GetLocG(vector<Node2D>& element, vector<vector<double>>& locG)
//{
//	for (int i = 0; i < 4; i++)
//	{
//		for (int j = 0; j < 4; j++)
//		{
//			locG[i][j] = 0;
//		}
//	}
//	for (int p_id = 0; p_id < 9; p_id++)
//	{
//		vector<vector<double>> J_1(2, vector<double>(2));
//		double det_J = CalcJ_1(J_1, p_id, element);
//
//		// intagrate_points_2D[p_id].derivative[i] - градиент i-ой базисной функции в точке p_id
//		vector<vector<double>> modifided_grad(4, vector<double>(2));
//		double modifided_coef = intagrate_points_2D[p_id].koef * abs(det_J);
//		for (int i = 0; i < 4; i++)
//		{
//			// modifided_grad[i] = J_1 * intagrate_points_2D[p_id].derivative[i]
//			for (int j = 0; j < 2; j++)
//			{
//				for (int k = 0; k < 2; k++)
//				{
//					modifided_grad[i][j] += J_1[j][k] * intagrate_points_2D[p_id].derivative[i][k];
//				}
//			}
//		}
//
//		for (int i = 0; i < 4; i++)
//		{
//			for (int j = 0; j < 4; j++)
//			{
//				locG[i][j] += mult(modifided_grad[i], modifided_grad[j]) * modifided_coef;
//			}
//		}
//
//	}
//
//	return 0;
//}
//
//int Integrate_Gauss3Method2D::GetLocM(vector<Node2D>& element, vector<vector<double>>& locM)
//{
//	for (int i = 0; i < 4; i++)
//	{
//		for (int j = 0; j < 4; j++)
//		{
//			locM[i][j] = 0;
//		}
//	}
//	for (int p_id = 0; p_id < 9; p_id++)
//	{
//		double det_J = CalcdetJ(p_id, element);
//
//		// intagrate_points_2D[p_id].basis_func[i] - i-а€ базисна€ функции в точке p_id
//		double modifided_coef = intagrate_points_2D[p_id].koef * abs(det_J);
//		for (int i = 0; i < 4; i++)
//		{
//			for (int j = 0; j < 4; j++)
//			{
//				locM[i][j] += intagrate_points_2D[p_id].basis_func[i] * intagrate_points_2D[p_id].basis_func[j] * modifided_coef;
//			}
//		}
//	}
//
//	return 0;
//}


