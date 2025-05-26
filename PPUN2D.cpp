#include "/usr/local/apps/tecplot2022r1/360ex_2022r1/include/TECIO.h"
#include <iostream>
#include <fstream> 
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <string>
#include <cstdlib>
#include <math.h>
#include <unordered_set>
#include <chrono>

using namespace std;

auto start_total = std::chrono::high_resolution_clock::now(); // Start to calculate total execution time
auto start_single = std::chrono::high_resolution_clock::now(); // Start to calculate execution time of one loop on Xw and previous calcuation 

int main(int argc, char** argv) {

    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " config.cas input_file.szplt output_file.szplt" << std::endl;
        return 1;
    }

    // Gas constants (ideal gas model)
    const double R_air = 287.058;
    const double gamma_air = 1.4;
    
    // Variables to be read in a configuration file (.cas)
    int momentumX_id = 4;
	int momentumY_id = 5;
    int Vx_id = 4;
    int Vy_id = 5;
	int rho_id = 3;
	int p_id = 9;
	int T_id = 10;
	int mut_id = 18;
	int k_id = 7;
    double AoA = 0;
    double rho_inf = 1.225;
	double T_inf = 288.15;
	double V_inf = 100;
	double chord = 1.0;
    double x_min = 0.0; // control volume
    double x_max = 0.0; // control volume 
    double y_min = 0.0; // control volume
    double y_max = 0.0; // control volume

    // Variables for region selection (drag-breakdown)
    double K_bl  = 1.0; // cut-off value viscous sensor
    double K_w   = 1.0; // cut-off value shock wave sensor
    double omega_tsh = 1.0; // cut off-value omega sensor
    int N_wall_margin = 0; // upper and lower margin for viscous region
    int N_Bl_margin = 0; // upper and lower margin for viscous BL region
    int N_wake_margin = 0; // upper and lower margin for viscous wake region
    int N_shock_margin = 0; // upper and lower margin for shock region
    int N_wv_bl_interaction = 0; // margin for shock wave - BL interface
    int N_shock_streamwise_margin = 0; // margin for shock region extension in streamwise direction
    double max_y_bl_udf = 0.01; // UDF max y height BL region 
    double min_y_bl_udf = -0.01; // UDF min y height BL region
    double wake_slope = 1.0;

    // Pole coordinates in Lamb vector-based methods
    double x_pole = 0.0;
    double y_pole = 0.0;

    // Asintotic flow field variables
    double p_inf = 0.0;
	double M_inf = 0.0;
    double q_inf = 0.0;

    // Flags
    bool split = false; // if split = 1 then quadrilateral cells are splitted in two triangular cells 
    bool divTauOff = false; // if divTauOff == 1 then div(tau) term is not considered in the computation of DV drag
    bool vrt_std_on = false; // if vrt_std_on then Drag and Lift are computed also with the standard Lamb-vector based method

    // Constant Near-field computed quantities
    double Cd_nf_stored = 0.0;
    double Cl_nf_stored = 0.0;

    // Read the configuration file
   	ifstream case_file(argv[2]);
	if (case_file.is_open()) {
		string line;
		int i = 0;
		int m = 0;
		while (getline(case_file, line)) {
			i = line.find_last_of("\t ");
			m = line.find_first_of("\t ");
			if (line.substr(0, m) == "AoA") {
				AoA = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "RHO_ID") {
				rho_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "RVX_ID") {
				momentumX_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "RVY_ID") {
				momentumY_id = atof(line.substr(i + 1).data());
			}
            else if (line.substr(0, m) == "VX_ID") {
			    Vx_id = atof(line.substr(i + 1).data());
			}
            else if (line.substr(0, m) == "VY_ID") {
				Vy_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "K_ID") {
				k_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "P_ID") {
				p_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "T_ID") {
				T_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "MUT_ID") {
				mut_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "RHO_INF") {
				rho_inf = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "T_INF") {
				T_inf = atof(line.substr(i + 1).data());
			}
            else if (line.substr(0, m) == "V_INF") {
				V_inf = atof(line.substr(i + 1).data());
			}
            else if (line.substr(0, m) == "AoA") {
				AoA = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "X_MIN") {
				x_min = atof(line.substr(i + 1).data());
			}
            else if (line.substr(0, m) == "X_MAX") {
				x_max = atof(line.substr(i + 1).data());
			}
            else if (line.substr(0, m) == "Y_MIN") {
				y_min = atof(line.substr(i + 1).data());
			}
            else if (line.substr(0, m) == "Y_MAX") {
				y_max = atof(line.substr(i + 1).data());
			}
            else if (line.substr(0, m) == "K_bl") {
				K_bl = atof(line.substr(i + 1).data());
			}
            else if (line.substr(0, m) == "K_w") {
				K_w = atof(line.substr(i + 1).data());
			}
            else if (line.substr(0, m) == "SPLIT") {
				split = atof(line.substr(i + 1).data());
			} 
            else if (line.substr(0, m) == "DIV_TAU_OFF") {
				divTauOff = atof(line.substr(i + 1).data());
			}
            else if (line.substr(0, m) == "VRT_STD_ON") {
				vrt_std_on = atof(line.substr(i + 1).data());   
			}
            else if (line.substr(0, m) == "OMEGA_TSH") {
				omega_tsh = atof(line.substr(i + 1).data());   
			} 
            else if (line.substr(0, m) == "N_WALL_MARGIN") {
				N_wall_margin = atof(line.substr(i + 1).data());   
			} 
            else if (line.substr(0, m) == "N_SHOCK_MARGIN") {
				N_shock_margin = atof(line.substr(i + 1).data());   
			}
            else if (line.substr(0, m) == "N_BL_MARGIN") {
				N_Bl_margin = atof(line.substr(i + 1).data());   
			} 
            else if (line.substr(0, m) == "N_WV_BL_INTERFACE") {
				N_wv_bl_interaction = atof(line.substr(i + 1).data());   
			}
            else if (line.substr(0, m) == "N_SHOCK_STREAMWISE_MARGIN") {
				N_shock_streamwise_margin = atof(line.substr(i + 1).data());   
			}
            else if (line.substr(0, m) == "N_WAKE_MARGIN") {
				N_wake_margin = atof(line.substr(i + 1).data());   
			} 
            else if (line.substr(0, m) == "WAKE_SLOPE") {
				wake_slope = atof(line.substr(i + 1).data());   
			} 
            else if (line.substr(0, m) == "MAX_Y_BL") {
				max_y_bl_udf = atof(line.substr(i + 1).data());   
			} 
            else if (line.substr(0, m) == "MIN_Y_BL") {
				min_y_bl_udf = atof(line.substr(i + 1).data());   
			} 
            else if (line.substr(0, m) == "X_POLE") {
				x_pole = atof(line.substr(i + 1).data());   
			}
            else if (line.substr(0, m) == "Y_POLE") {
				y_pole = atof(line.substr(i + 1).data());   
			}
		}
		case_file.close();
	}

    // Asymptotic flow field
    p_inf= rho_inf * R_air * T_inf;
	M_inf = V_inf / (sqrt(gamma_air * R_air * T_inf));
    q_inf = 0.5 * rho_inf * pow(V_inf,2.0);
    std::cout << "Mach_inf = " << M_inf << endl;
    
    // **************************************** READ FILE (.szplt) ***********************************************
  
    // Open input file
    void* inputFileHandle = NULL;
    int32_t R = tecFileReaderOpen(argv[1], &inputFileHandle);
    if ( R != 0) throw std::runtime_error("Failed to open input file.");

    // General dataset informations
    char* dataSetTitle = NULL;
    R = tecDataSetGetTitle(inputFileHandle, &dataSetTitle);
    int32_t numVars = 0;   // original .szplt number of variables
    int32_t nodeVars = 6;  // new nodal variables

    // Enumeration to link each variable to the correct index 
    // Cell centered variables
    enum CellVarsIndex {
        CENTROID_X, CENTROID_Y, CELL_AREA, NORMAL_X1, NORMAL_Y1, NORMAL_X2, NORMAL_Y2, NORMAL_X3, NORMAL_Y3, NORMAL_X4, NORMAL_Y4, DU_DX, DU_DY, DV_DX, DV_DY,
        TAU11_NF, TAU12_NF, TAU21_NF, TAU22_NF, TAU11_FF, TAU12_FF, TAU21_FF, TAU22_FF, CELL_TYPE,
        CELL_FLUX_DRAG_PT, CELL_FLUX_DRAG_DV, CELL_FLUX_DRAG_PT_WV, CELL_FLUX_DRAG_PT_VS, CELL_FLUX_DRAG_PT_SP, CELL_FLUX_DRAG_DV_WV, CELL_FLUX_DRAG_DV_VS, CELL_FLUX_DRAG_DV_SP,
        DRHO_DX, DRHO_DY, DS_DX, DS_DY, DH_DX, DH_DY, DK_DX, DK_DY, DKT_DX, DKT_DY, DDY_RHODKTDX, DDX_RHODKTDY, 
        FS_STD_WV, FS_THD_WV, FS_STD_VS, FS_THD_VS, FS_STD_SP, FS_THD_SP, FS_STD, FS_THD, LAMB_STD_X, LAMB_STD_Y, LAMB_THD_X, LAMB_THD_Y,
        NUM_CELL_VARS
    };
    // Cell centered temporary variables
    enum CellTempVarsIndex {
        DTAU11DX_NF, DTAU12DY_NF, DTAU12DX_NF, DTAU22DY_NF, DTAU11DX_FF, DTAU12DY_FF, DTAU12DX_FF, DTAU22DY_FF, DP_DX, DP_DY, U_CC, V_CC, T_CC, RHO_CC, MU_D_CC, MU_T_CC, K_CC, R_XC, R_YC,
        R_DOT_L_std, R_DOT_L_thd, R_RHOL_std_XX, R_RHOL_std_XY, R_RHOL_std_YX, R_RHOL_std_YY, 
        R_RHOL_thd_XX, R_RHOL_thd_XY, R_RHOL_thd_YX, R_RHOL_thd_YY, GRAD_RDOTL_std_X, GRAD_RDOTL_std_Y, 
        GRAD_RDOTL_thd_X, GRAD_RDOTL_thd_Y, DIV_R_RHOL_std_X, DIV_R_RHOL_std_Y, DIV_R_RHOL_thd_X, DIV_R_RHOL_thd_Y,
        NUM_TEMP_CELL_VARS
    };

    int32_t cellVars = NUM_CELL_VARS; // new cell centered variables
    int32_t cellTempVars = NUM_TEMP_CELL_VARS; // cell centered temp variables --> 13 variabili ma per ora la posizione cellTempVars - 1 non è occcupata!! => 12 attive
    R = tecDataSetGetNumVars(inputFileHandle, &numVars);
    int32_t numVarsOut = numVars + nodeVars + cellVars;
    //int32_t numVarsOut = numVars + nodeVars + cellVars + cellTempVars; // total output variables 

    // Variables list 
    std::ostringstream varListStream;
    for (int32_t var = 1; var <= numVars; ++var) {
        char* varName = NULL;
        R = tecVarGetName(inputFileHandle, var, &varName);
        varListStream << varName;
        if (var < numVars) varListStream << ',';
        tecStringFree(&varName);
    }

    // Add new variables (same order of 'enum' !)
    varListStream << ",Entropy_variation,Entalpy_variation,nodeID_ext,nodeID_body,r_x,r_y,centroid_x,centroid_y,cell_area,nx1,ny1,nx2,ny2,nx3,ny3,nx4,ny4,du_dx,du_dy,dv_dx,dv_dy,tau11_nf,tau12_nf,tau21_nf,tau22_nf,tau11_ff,tau12_ff,tau21_ff,tau22_ff,cell_type,flux_drag_PT,flux_drag_DV,flux_drag_PT_wave,flux_drag_PT_viscous,flux_drag_PT_spurious,flux_drag_DV_wave,flux_drag_DV_viscous,flux_drag_DV_spurious,DrhoDx,DrhoDy,dSdx,dSdy,dHdx,dHdy,dKdx,dKdy,dKtdx,dKtdy,DDY_RHODKTDX,DDX_RHODKTDY,FS_STD_WV,FS_THD_WV,FS_STD_VS,FS_THD_VS,FS_STD_SP,FS_THD_SP,FS_STD,FS_THD,LAMB_STD_X, LAMB_STD_Y, LAMB_THD_X, LAMB_THD_Y";
    //varListStream << ",Entropy_variation,Entalpy_variation,nodeID_ext,nodeID_body,r_x,r_y,centroid_x,centroid_y,cell_area,nx1,ny1,nx2,ny2,nx3,ny3,nx4,ny4,du_dx,du_dy,dv_dx,dv_dy,tau11_nf,tau12_nf,tau21_nf,tau22_nf,tau11_ff,tau12_ff,tau21_ff,tau22_ff,cell_type,flux_drag_PT,flux_drag_DV,flux_drag_PT_wave,flux_drag_PT_viscous,flux_drag_PT_spurious,flux_drag_DV_wave,flux_drag_DV_viscous,flux_drag_DV_spurious,DrhoDx,DrhoDy,dSdx,dSdy,dHdx,dHdy,dKdx,dKdy,dKtdx,dKtdy,DDY_RHODKTDX,DDX_RHODKTDY,dTau11dx_nf,dTau12dy_nf,dTau12dx_nf,dTau22dy_nf,dTau11dx_ff,dTau12dy_ff,dTau12dx_ff,dTau22dy_ff,dPdX,dPdY,u_cc,v_cc,T_cc,rho_cc,mu_d_cc,mu_t_cc,k_cc,r_xc,r_yc";
    
    int32_t numZones = 0;
    R = tecDataSetGetNumZones(inputFileHandle, &numZones);
    if (numZones <= 0) {
        throw std::runtime_error("Invalid number of zones: " + std::to_string(numZones));
    }

    // Adapt types to TecIO library
     int32_t isDouble = 0;
    int32_t const FieldDataType_Double = 2;
    if (numZones > 0)
    {
        int32_t type = 0;
        R = tecZoneVarGetType(inputFileHandle, 1, 1, &type);
        if (type == FieldDataType_Double)
            isDouble = 1;
    } 
    
    // Zone info
    char* zoneTitle = NULL;
    R = tecZoneGetTitle(inputFileHandle, 1, &zoneTitle);

    int32_t zoneType = 0;
    R = tecZoneGetType(inputFileHandle, 1, &zoneType);
    if (zoneType == 6 || zoneType == 7)
        throw std::runtime_error("Unsupported zone type.");

    int64_t iMax = 0, jMax = 0, kMax = 0;
    tecZoneGetIJK(inputFileHandle, 1, &iMax, &jMax, &kMax);

    // For unstructured grids
    int64_t numNodes = iMax;     // Number of nodes
    int64_t numCells = jMax;     // Number of cells (elements)
    int64_t nodesPerCell = 4;    // Each cell has 4 nodes; for triangular cells one node is repeated!

    std::cout << "Zone info:\n" << "Number of Nodes: " << numNodes << ", Number of Cells: " << numCells << std::endl;

    // Variables informations
    std::vector<int32_t> varTypes(numVars);
    std::vector<int32_t> passiveVarList(numVars);
    std::vector<int32_t> valueLocation(numVars);
    std::vector<int32_t> shareVarFromZone(numVars);
    
    for (int32_t var = 1; var <= numVars; ++var)
    {
        R = tecZoneVarGetType(inputFileHandle, 1, var, &varTypes[var - 1]);
        R = tecZoneVarGetSharedZone(inputFileHandle, 1, var, &shareVarFromZone[var - 1]);
        R = tecZoneVarGetValueLocation(inputFileHandle, 1, var, &valueLocation[var - 1]);
        R = tecZoneVarIsPassive(inputFileHandle, 1, var, &passiveVarList[var - 1]);
    }

    int32_t shareConnectivityFromZone;
    R = tecZoneConnectivityGetSharedZone(inputFileHandle, 1, &shareConnectivityFromZone);

    int32_t faceNeighborMode;
    R = tecZoneFaceNbrGetMode(inputFileHandle, 1, &faceNeighborMode);

    int64_t numFaceConnections;
    R = tecZoneFaceNbrGetNumConnections(inputFileHandle, 1, &numFaceConnections);

    // Creation of data structures to manipulate data
    std::vector < std::vector < double> > nodeVariables(numNodes, std::vector<double>(numVars + nodeVars , 0.0));         // Nodal variables
    std::vector < std::vector < double> > cellCenteredVariables(numCells, std::vector<double>(NUM_CELL_VARS, 0.0));       // Cell centered variables
    std::vector < std::vector < double> > ccVariables(numCells, std::vector<double>(NUM_TEMP_CELL_VARS, 0.0));            // Cell centered temp variables
    std::vector < std::vector <int32_t> > connectivity(numCells, std::vector<int32_t>(nodesPerCell, 0));                  // Connectivity list (32 bit)    
    vector<vector<double>> normals(numCells, vector<double>(nodesPerCell * 2, 0.0));                                      // Normals matrix (NumCells x 8)

    // Reading node variables
    for (int32_t var = 1; var <= numVars; var++) {
            int64_t numValues = 0;
            R = tecZoneVarGetNumValues(inputFileHandle, 1, var, &numValues);
            if (numValues != numNodes) {
                throw std::runtime_error("Mismatch in variable size.");
            }
            // Lettura dei valori in un vettore temporaneo (formato DATAPACKING=BLOCK)
            std::vector<double> tempValues(numNodes, 0.0); // vettore temporaneo per i valori delle variabili
            R = tecZoneVarGetDoubleValues(inputFileHandle, 1, var, 1, numNodes, &tempValues[0]);
            // Copia dei valori nella matrice nodeVariables
            for(int64_t node = 0; node < numNodes; node++) {
                nodeVariables[node][var - 1] = tempValues[node]; // riempi ogni colonna corrispondente alla variabile 'var' col vettore 'tempValues' letto dal file .szplt
            }
    }

    // Reading adjacent faces, if any
    if (numFaceConnections > 0)
    {
        int64_t numFaceValues;
        R = tecZoneFaceNbrGetNumValues(inputFileHandle, 1, &numFaceValues);
        int32_t are64Bit;
        R = tecZoneFaceNbrsAre64Bit(inputFileHandle, 1, &are64Bit);
        if (are64Bit)
        {
            std::vector<int64_t> faceConnections(numFaceValues);
            R = tecZoneFaceNbrGetConnections64(inputFileHandle, 1, &faceConnections[0]);
        }
        else
        {
            std::vector<int32_t> faceConnections(numFaceValues);
            R = tecZoneFaceNbrGetConnections(inputFileHandle, 1, &faceConnections[0]);
        }
    }

    // Reading connectivity list
    if (zoneType != 0 && shareConnectivityFromZone == 0)
    {
        int64_t numValues = 0;
        R = tecZoneNodeMapGetNumValues(inputFileHandle, 1, numCells, &numValues);
        std::vector<int32_t> nodeMap(numValues);
        R = tecZoneNodeMapGet(inputFileHandle, 1, 1, numCells, &nodeMap[0]);
        for (int64_t elem = 0; elem < numCells; elem++) {
            for(int localNode = 0; localNode < nodesPerCell; localNode++) {
                connectivity[elem][localNode] = nodeMap[elem * nodesPerCell + localNode];
            }
        }                  
    }

    std::cout << "NodeVariables matrix : " << nodeVariables.size() << " x " << nodeVariables[0].size() << std::endl;
    std::cout << "CellCenteredVariables matrix : " << cellCenteredVariables.size() << " x " << cellCenteredVariables[0].size() << std::endl;
    std::cout << "Connectivity matrix: " << connectivity.size() << " x " << connectivity[0].size() << "\n" << std::endl;

    // ******************************** DATA MANIPULATION (POST-PROCESSING) ******************************************

    // +++++++++++++++++++++++++ PART 1: QUANTITIES UNDEPENDENT ON THE CONTROL VOLUME ++++++++++++++++++++++++

    // Split quadrilateral cell in two triangular cells => NewNumCells = N_triang + 2xN_quad
    // Count quadrilateral and triangular cells
    int numQuadCells = 0;
    int numTriangCells = 0;
    for (int cell = 0; cell < numCells; cell++) {
        // Extract nodes from connectivity matrix
        int32_t n1 = connectivity[cell][0];
        int32_t n2 = connectivity[cell][1];
        int32_t n3 = connectivity[cell][2];
        int32_t n4 = connectivity[cell][3];

        if (n3 == n4) {
        numTriangCells++;
        } else {
        numQuadCells++;
        }
    }

    // If split == 1 then split quadrilateral cells
    if (split) {

        // Calculate the new number of cells and the new connectivity matrix
        int newNumCells = (numQuadCells * 2) + numTriangCells;
        vector<vector<int32_t>> newConnectivity;

        // Iterate over the original cells and add the first triangle
        for (int cell = 0; cell < numCells; cell++) {

            int32_t n1 = connectivity[cell][0];
            int32_t n2 = connectivity[cell][1];
            int32_t n3 = connectivity[cell][2];
            int32_t n4 = connectivity[cell][3];

            if (n3 == n4) { // triangular => copy

                newConnectivity.push_back({n1,n2,n3,n3});

            } else { // quadrilateral => split

                newConnectivity.push_back({n1,n2,n3,n3});
                newConnectivity.push_back({n3,n4,n1,n1});
            } 
        }

        // Replace the old data structure with the new one
        connectivity = std::move(newConnectivity);
        numCells = newNumCells;
        cellCenteredVariables.resize(numCells, vector<double>(cellVars, 0.0));
        ccVariables.resize(numCells, vector<double>(cellTempVars, 0.0));
        normals.resize(numCells, vector<double>(nodesPerCell * 2, 0.0));


        std::cout << "Number of original quadrilateral cells: " << numQuadCells << std::endl;
        std::cout << "Number of original triangular cells: " << numTriangCells << std::endl;
        std::cout << "Number of cells after splitting: " << numCells << std::endl;
        std::cout << "New connectivity matrix: " << connectivity.size() << " x " << connectivity[0].size() << "\n" << std::endl; 
        std::cout << "New cellCentered matrix: " << cellCenteredVariables.size() << " x " << cellCenteredVariables[0].size() << "\n" << std::endl; 
        std::cout << "New normals matrix: " << normals.size() << " x " << normals[0].size() << "\n" << std::endl; 
    }

    // Rotate the grid when AoA != 0 in order to align x and y axes to wind axes. In this way aerodynamic forces will be computed in an easier way.
    double AoA_rad = AoA * M_PI / 180.0; // converts alpha from deg to rad

    if (AoA!=0) {	// To save computational time on zero AoA cases.
        // Counterclockwise rotation
        double cos_alpha = cos(-AoA_rad); 
        double sin_alpha = sin(-AoA_rad); 

        for (int node = 0; node < numNodes; node++) {
            // Rotate coordinates
            double x = nodeVariables[node][0];
            double y = nodeVariables[node][1];

            double x_rotated = x * cos_alpha - y * sin_alpha;
            double y_rotated = x * sin_alpha + y * cos_alpha;

            nodeVariables[node][0] = x_rotated;
            nodeVariables[node][1] = y_rotated;

            // Rotate momentum components
            double momentumX = nodeVariables[node][momentumX_id - 1];
            double momentumY = nodeVariables[node][momentumY_id - 1];

            double momentumX_rotated = momentumX * cos_alpha - momentumY * sin_alpha;
            double momentumY_rotated = momentumX * sin_alpha + momentumY * cos_alpha;

            nodeVariables[node][momentumX_id - 1] = momentumX_rotated;
            nodeVariables[node][momentumY_id - 1] = momentumY_rotated;
        }
    }

    // Centroids and cell areas
    for(int64_t cell = 0; cell < numCells; cell++) {
        int numFaces = 0;
        if(connectivity[cell][2] == connectivity[cell][3]) {
            numFaces = 3;
        }
        else{
            numFaces = 4;
        }

        double xSum = 0.0, ySum = 0.0;
        for (int localNode = 0; localNode < numFaces; localNode++) {
            int64_t node = connectivity[cell][localNode] - 1; // tecIO is 1-based, shift
            xSum += nodeVariables[node][0];
            ySum += nodeVariables[node][1];
        }

        // Centroids
        cellCenteredVariables[cell][CENTROID_X] = xSum / numFaces; // centroid_x
        cellCenteredVariables[cell][CENTROID_Y] = ySum / numFaces; // centroid_y

        // Areas
        double area = 0.0;
        if (numFaces == 3) { // triangular face --> cartesian coordinates area formula A = 1/2*|x1(y2-y3) + x2(y3-y1) + x3(y1-y2)|
            int64_t n1 = connectivity[cell][0] - 1;
            int64_t n2 = connectivity[cell][1] - 1;
            int64_t n3 = connectivity[cell][2] - 1;
            area = 0.5 * fabs(
                nodeVariables[n1][0] * (nodeVariables[n2][1] - nodeVariables[n3][1]) +
                nodeVariables[n2][0] * (nodeVariables[n3][1] - nodeVariables[n1][1]) +
                nodeVariables[n3][0] * (nodeVariables[n1][1] - nodeVariables[n2][1])
            );
        }
        else { // quadrilateral faces --> same formula but considering each rectangle splitted in two triangle
            int64_t n1 = connectivity[cell][0] - 1;
            int64_t n2 = connectivity[cell][1] - 1;
            int64_t n3 = connectivity[cell][2] - 1;
            int64_t n4 = connectivity[cell][3] - 1; 
            double area1 = 0.5 * fabs(
                nodeVariables[n1][0] * (nodeVariables[n2][1] - nodeVariables[n3][1]) +
                nodeVariables[n2][0] * (nodeVariables[n3][1] - nodeVariables[n1][1]) +
                nodeVariables[n3][0] * (nodeVariables[n1][1] - nodeVariables[n2][1])
            );
            double area2 = 0.5 * fabs(
                nodeVariables[n1][0] * (nodeVariables[n3][1] - nodeVariables[n4][1]) +
                nodeVariables[n3][0] * (nodeVariables[n4][1] - nodeVariables[n1][1]) +
                nodeVariables[n4][0] * (nodeVariables[n1][1] - nodeVariables[n3][1])
            );
            area = area1 + area2;
        }
        cellCenteredVariables[cell][CELL_AREA] = area; 
    }

    /* Face centered non unitary normals (n*dS) */
    // Loop on cells and faces
    for (int64_t cell = 0; cell < numCells; cell++)  {
        int maxNumFaces = 4; 
        double n_x = 0.0; // normal vector - x component
        double n_y = 0.0; // normal vector - y component
        double orth = 0;  // cell orthogonal vector = t1 x t2
        int64_t node1, node2, node3;
        for (int localFace = 0; localFace < maxNumFaces; localFace++) {
            // Compute orth vector
            node1 = connectivity[cell][0];
            node2 = connectivity[cell][1];
            node3 = connectivity[cell][2];
            double x1 = nodeVariables[node1 - 1][0];
            double y1 = nodeVariables[node1 - 1][1];
            double x2 = nodeVariables[node2 - 1][0];
            double y2 = nodeVariables[node2 - 1][1];
            double x3 = nodeVariables[node3 - 1][0];
            double y3 = nodeVariables[node3 - 1][1];
            double t1_x = x2 - x1;
            double t1_y = y2 - y1;
            double t2_x = x3 - x2;
            double t2_y = y3 - y2;
            orth = t1_x*t2_y-t1_y*t2_x;

            if(localFace == 0) { // face 1
                node1 = connectivity[cell][0];  
                node2 = connectivity[cell][1];  
            } else if (localFace == 1) { // face 2
                node1 = connectivity[cell][1]; 
                node2 = connectivity[cell][2]; 
            } else if (localFace == 2) { // face 3
                node1 = connectivity[cell][2]; 
                node2 = connectivity[cell][3]; 
            } else if (localFace == 3) { // face 4
                node1 = connectivity[cell][3]; 
                node2 = connectivity[cell][0]; 
            }
            
            // Nodes coordinates (reminder: tecIO is 1-based!)
            x1 = nodeVariables[node1 - 1][0];
            y1 = nodeVariables[node1 - 1][1];
            x2 = nodeVariables[node2 - 1][0];
            y2 = nodeVariables[node2 - 1][1];

            // Face tangent vector
            double dl_x = x2 - x1;
            double dl_y = y2 - y1;

            // Compute normal vector components 
            double n_x = 0;
            double n_y = 0;
            if (orth > 0) { // anti-clockwise sorting of cell nodes 
                n_x = dl_y;
                n_y = - dl_x;
            } else { // clockwise sorting of cell nodes
                n_x = - dl_y;
                n_y = dl_x;
            } 

            // Save normal components in 'normals'
            normals[cell][localFace * 2] = n_x; // n_x
            normals[cell][localFace * 2 + 1] = n_y; // n_y
        }
    }

    /* Compute faces shared by adjacent cells */
    // Map: key = face defined by two nodes; value = list of cells (index) wich includes that face.
    map<pair<int32_t, int32_t>, vector<int>> faceToCells; // faceToCells associates one face to a list of cells which include that face
    // Map population
    for (int cell = 0; cell < numCells; cell++){
        int32_t node1 = connectivity[cell][0]; 
        int32_t node2 = connectivity[cell][1];
        int32_t node3 = connectivity[cell][2];
        int32_t node4 = connectivity[cell][3];

        bool isTriangle = (node3 == node4);

        vector <pair<int32_t, int32_t>> faces;
        if (isTriangle) {
            faces.push_back({min(node1, node2), max(node1, node2)});
            faces.push_back({min(node2, node3), max(node2, node3)});
            faces.push_back({min(node3, node1), max(node3, node1)});
        } else {
            faces.push_back({min(node1, node2), max(node1, node2)});
            faces.push_back({min(node2, node3), max(node2, node3)});
            faces.push_back({min(node3, node4), max(node3, node4)});
            faces.push_back({min(node4, node1), max(node4, node1)});
        }

        for (const auto& face : faces) {
            faceToCells[face].push_back(cell);
        }
    }
    
    // Compute cell centered velocity and pressure gradients (GREEN-GAUSS nodes based method)
    for(int cell = 0; cell < numCells; cell++) {

        double A_cell = cellCenteredVariables[cell][CELL_AREA]; // cell area
        double grad_u[2] = {0.0, 0.0}; // du/dx, du/dy
        double grad_v[2] = {0.0, 0.0}; // dv/dx, dv/dy
        double grad_p[2] = {0.0, 0.0}; // dp/dx, dp/dy

        // Compute gradient contributions on faces
        for(int localFace = 0; localFace < 4; localFace++) {

            // Face nodes index
            int64_t node1 = connectivity[cell][localFace];
            int64_t node2 = connectivity[cell][(localFace + 1) % 4];
            
            if(node1 == node2) { // triangular cell => 3 faces
                continue;
            }

            // Nodal values
            double u1    = nodeVariables[node1 - 1][momentumX_id - 1] / nodeVariables[node1 - 1][rho_id - 1];
            double u2    = nodeVariables[node2 - 1][momentumX_id - 1] / nodeVariables[node2 - 1][rho_id - 1];
            double v1    = nodeVariables[node1 - 1][momentumY_id - 1] / nodeVariables[node1 - 1][rho_id - 1];
            double v2    = nodeVariables[node2 - 1][momentumY_id - 1] / nodeVariables[node2 - 1][rho_id - 1];
            double p1    = nodeVariables[node1 - 1][p_id - 1];
            double p2    = nodeVariables[node2 - 1][p_id - 1];
            
            // Face centered interpolation
            double u_face = 0.5 * (u1 + u2);
            double v_face = 0.5 * (v1 + v2);
            double p_face = 0.5 * (p1 + p2);

            // Normals
            double n_x = normals[cell][localFace * 2];
            double n_y = normals[cell][localFace * 2 + 1];

            // Sum on faces
            grad_u[0] += u_face * n_x;
            grad_u[1] += u_face * n_y;
            grad_v[0] += v_face * n_x;
            grad_v[1] += v_face * n_y;
            grad_p[0] += p_face * n_x;
            grad_p[1] += p_face * n_y;
        }

        // Cell centered gradient components
        grad_u[0] = grad_u[0] / A_cell;
        grad_u[1] = grad_u[1] / A_cell;
        grad_v[0] = grad_v[0] / A_cell;
        grad_v[1] = grad_v[1] / A_cell;
        grad_p[0] = grad_p[0] / A_cell;
        grad_p[1] = grad_p[1] / A_cell;

        // Saving gradient components 
        cellCenteredVariables[cell][DU_DX]  = grad_u[0];  // du/dx
        cellCenteredVariables[cell][DU_DY]  = grad_u[1];  // du/dy 
        cellCenteredVariables[cell][DV_DX]  = grad_v[0];  // dv/dx
        cellCenteredVariables[cell][DV_DY]  = grad_v[1];  // dv/dy
        ccVariables[cell][DP_DX]            = grad_p[0];  // dp/dx
        ccVariables[cell][DP_DY]            = grad_p[1];  // dp/dy
    }

    /* Compute auxiliary cell centered values */
    for (int cell = 0; cell < numCells; cell++) {
        int numFaces = 0.0;
        double u_cc = 0.0, v_cc = 0.0, T_cc = 0.0, rho_cc = 0.0, mu_d_cc = 0.0, mu_t_cc = 0.0, k_cc = 0.0;
        for (int node = 0; node < 4; node++) {

                int64_t node1 = connectivity[cell][0];
                int64_t node2 = connectivity[cell][1];
                int64_t node3 = connectivity[cell][2];
                int64_t node4 = connectivity[cell][3];

                if (node3 == node4) {
                    numFaces = 3;
                } else {
                    numFaces = 4;
                }

                double T1    = nodeVariables[node1 - 1][T_id - 1];
                double T2    = nodeVariables[node2 - 1][T_id - 1];
                double T3    = nodeVariables[node3 - 1][T_id - 1];
                double T4 = 0;
                if (numFaces==4) {
                    T4 = nodeVariables[node4 - 1][T_id - 1];
                }

                double mu_d1 = 1.458e-6 * sqrt(T1) / (1 + 110.4 / T1);
                double mu_d2 = 1.458e-6 * sqrt(T2) / (1 + 110.4 / T2);
                double mu_d3 = 1.458e-6 * sqrt(T3) / (1 + 110.4 / T3);
                double mu_d4 = 0;
                if (numFaces==4) {
                    mu_d4 = 1.458e-6 * sqrt(T4) / (1 + 110.4 / T4);
                }

                double mu_t1 = nodeVariables[node1 - 1][mut_id - 1];
                double mu_t2 = nodeVariables[node2 - 1][mut_id - 1];
                double mu_t3 = nodeVariables[node3 - 1][mut_id - 1];
                double mu_t4 = 0;
                if (numFaces==4) {
                    mu_t4 = nodeVariables[node4 - 1][mut_id - 1];
                }

                double u1    = nodeVariables[node1 - 1][momentumX_id - 1] / nodeVariables[node1 - 1][rho_id - 1];
                double u2    = nodeVariables[node2 - 1][momentumX_id - 1] / nodeVariables[node2 - 1][rho_id - 1];
                double u3    = nodeVariables[node3 - 1][momentumX_id - 1] / nodeVariables[node3- 1][rho_id - 1];
                double u4 = 0;
                if (numFaces==4) {
                    u4    = nodeVariables[node4 - 1][momentumX_id - 1] / nodeVariables[node4 - 1][rho_id - 1];
                }

                double v1    = nodeVariables[node1 - 1][momentumY_id - 1] / nodeVariables[node1 - 1][rho_id - 1];
                double v2    = nodeVariables[node2 - 1][momentumY_id - 1] / nodeVariables[node2 - 1][rho_id - 1];
                double v3    = nodeVariables[node3- 1][momentumY_id - 1] / nodeVariables[node3 - 1][rho_id - 1];
                double v4 = 0;
                if (numFaces==4) {
                    v4    = nodeVariables[node4 - 1][momentumY_id - 1] / nodeVariables[node4 - 1][rho_id - 1];
                }

                double rho1  = nodeVariables[node1 - 1][rho_id - 1];
                double rho2  = nodeVariables[node2 - 1][rho_id - 1];
                double rho3  = nodeVariables[node3 - 1][rho_id - 1];
                double rho4 = 0;
                if (numFaces==4) {
                    rho4  = nodeVariables[node4 - 1][rho_id - 1];
                }

                double k1    = nodeVariables[node1 - 1][k_id - 1];
                double k2    = nodeVariables[node2 - 1][k_id - 1];
                double k3    = nodeVariables[node3 - 1][k_id - 1];
                double k4 = 0;
                if (numFaces==4) {
                    k4    = nodeVariables[node4 - 1][k_id - 1];
                }
                
                // Sum node values
                T_cc = T1 + T2 + T3 + T4;
                mu_d_cc = mu_d1 + mu_d2 + mu_d3 + mu_d4;
                mu_t_cc = mu_t1 + mu_t2 + mu_t3 + mu_t4;
                u_cc = u1 + u2 + u3 + u4;
                v_cc = v1 + v2 + v3 + v4;
                rho_cc = rho1 + rho2 + rho3 + rho4;
                k_cc = k1 + k2 + k3 + k4;
            }

            // Divide for NumFaces
            //std:: cout << "Number of faces: " << numFaces << endl;
            T_cc /= numFaces;
            mu_d_cc /= numFaces;
            mu_t_cc /= numFaces;
            u_cc /= numFaces;
            v_cc /= numFaces;
            rho_cc /= numFaces;
            k_cc /= numFaces;

            // Save values in auxiliary structures
            ccVariables[cell][U_CC] = u_cc;
            ccVariables[cell][V_CC] = v_cc;
            ccVariables[cell][T_CC] = T_cc;
            ccVariables[cell][RHO_CC] = rho_cc;
            ccVariables[cell][MU_D_CC] = mu_d_cc;
            ccVariables[cell][MU_T_CC] = mu_t_cc;
            ccVariables[cell][K_CC] = k_cc;
    }

    // Compute viscous stress tensor: tau = 2*(mu_d + mu_t)(grad(V)_s0)
    for (int cell = 0; cell < numCells; cell++) {

        // Velocity gradient components
        double du_dx = cellCenteredVariables[cell][DU_DX];
        double du_dy = cellCenteredVariables[cell][DU_DY];
        double dv_dx = cellCenteredVariables[cell][DV_DX];
        double dv_dy = cellCenteredVariables[cell][DV_DY];

        double mu_d_cc = ccVariables[cell][MU_D_CC];
        double mu_t_cc = ccVariables[cell][MU_T_CC];
        double mu_cc = mu_d_cc + mu_t_cc;
    

        // Cell centered velocity field divergence 
        double divV = du_dx + dv_dy; // du/dx + dv/dy

        // Viscous stress tensor components
        double tau11  = 2 * mu_cc * (du_dx - (divV / 3)); // 2 * mu * (du/dx - divV/3)
        double tau12  = mu_cc * (dv_dx + du_dy);      // 2 * mu * ((dv/dx + du/dy)/2)
        //  double tau21 = tau12; (simmetry)
        double tau22  = 2 * mu_cc * (dv_dy - (divV / 3));  // 2 * mu * (dv/dy - divV/3) 

        // Saving tau components in cellCenteredVariables matrix
        cellCenteredVariables[cell][TAU11_NF] = tau11;
        cellCenteredVariables[cell][TAU12_NF] = tau12; 
        cellCenteredVariables[cell][TAU21_NF] = tau12;
        cellCenteredVariables[cell][TAU22_NF] = tau22;
    }

    // Compute components of tau gradient with Green-Gauss CELL BASED (tau_f = a*tau_p + (1-a)*tau_n)
    for (int cell = 0; cell < numCells; cell++) {

        double dTaudx[2] = {0.0, 0.0};
        double dTaudy[2] = {0.0, 0.0};
        double cellArea = cellCenteredVariables[cell][CELL_AREA];

        for (int localFace = 0; localFace < 4; localFace++) {
            int64_t node1 = connectivity[cell][localFace];
            int64_t node2 = connectivity[cell][(localFace + 1) % 4];

            if (node1 == node2) continue;

            double n_x = normals[cell][localFace * 2];
            double n_y = normals[cell][localFace * 2 + 1];

            // Face centroid
            double x_f = 0.5 * (nodeVariables[node1 - 1][0] + nodeVariables[node2 - 1][0]);
            double y_f = 0.5 * (nodeVariables[node1 - 1][1] + nodeVariables[node2 - 1][1]);

            // Looking for adjacent cell
            pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
            int adjacentCell = -1; // none
            auto it = faceToCells.find(face); // reference to faceToCells
            if (it != faceToCells.end() && it->second.size() == 2) { // check if 'face' is shared exactly by two cells
                adjacentCell = (it->second[0] == cell) ? it->second[1] : it->second[0]; // define adjacent cell (different from current cell)
                //cout << "FOUND N1: " << face.first << ", N2: " << face.second << ", WORKING CELL: " << cell << ", ADJACENT CELL: " << adjacentCell << "\n";
            }

            // Face centered tau interpolation
            double tau11_face = 0.0, tau12_face = 0.0, tau22_face = 0.0;

            if(adjacentCell != -1) { // if adjacent cell exists
                double tau11_n = cellCenteredVariables[cell][TAU11_NF];
                double tau11_p = cellCenteredVariables[adjacentCell][TAU11_NF];
                double tau12_n = cellCenteredVariables[cell][TAU12_NF];
                double tau12_p = cellCenteredVariables[adjacentCell][TAU12_NF];
                double tau22_n = cellCenteredVariables[cell][TAU22_NF];
                double tau22_p = cellCenteredVariables[adjacentCell][TAU22_NF];

                // Compute a 
                double a = 0.0;
                double dist_n = sqrt(pow(x_f - cellCenteredVariables[cell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[cell][CENTROID_Y], 2));
                double dist_p = sqrt(pow(x_f - cellCenteredVariables[adjacentCell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[adjacentCell][CENTROID_Y], 2));
                a = dist_n / (dist_p + dist_n);
                
                tau11_face = a * tau11_p + (1.0 - a) * tau11_n;
                tau12_face = a * tau12_p + (1.0 - a) * tau12_n;
                tau22_face = a * tau22_p + (1.0 - a) * tau22_n;
            } else {
                tau11_face = cellCenteredVariables[cell][TAU11_NF];
                tau12_face = cellCenteredVariables[cell][TAU12_NF];
                tau22_face = cellCenteredVariables[cell][TAU22_NF];
            }

            // Tau gradient components
            dTaudx[0] += tau11_face * n_x; // dtau11/dx
            dTaudy[1] += tau12_face * n_y; // dtau12/dy
            dTaudx[1] += tau12_face * n_x; // dtau12/dx
            dTaudy[0] += tau22_face * n_y; // dtau22/dy
        }

        dTaudx[0] /= cellArea; // dtau11/dx
        dTaudy[1] /= cellArea; // dtau12/dy
        dTaudx[1] /= cellArea; // dtau12/dx
        dTaudy[0] /= cellArea; // dtau22/dy
        
        ccVariables[cell][DTAU11DX_NF] = dTaudx[0];
        ccVariables[cell][DTAU12DY_NF] = dTaudy[1];
        ccVariables[cell][DTAU12DX_NF] = dTaudx[1];
        ccVariables[cell][DTAU22DY_NF] = dTaudy[0];
    }

    /* Quantities to compute drag in Thd Far-field methods: D = V_inf*INT_Omega[div(rho*f*V)dOmega] */

    // Coefficients fs1 e fs2
    double fs1 = 0, fs2 = 0;    
    fs1 = - 1.0 / (gamma_air*(pow(M_inf,2.0)));
    fs2 = - ((1.0 + (gamma_air - 1)*pow(M_inf,2.0)) / (2.0*pow(gamma_air,2.0)*pow(M_inf,4.0)));

    // Non-dimensional entropy variation on nodes (s-s_inf / R = ln[(rho_inf/rho)(T/T_inf)^(1/(gamma-1))])
    vector<double> deltaS_R(numNodes, 0.0); 
    for (int node = 0; node < numNodes; node++) {
        double rho = nodeVariables[node][rho_id - 1];
        double T   = nodeVariables[node][T_id - 1]; 
        double rho_ratio = rho_inf / rho;
        double T_ratio   = T / T_inf;
        deltaS_R[node] = log(rho_ratio * pow(T_ratio, 1.0 / (gamma_air - 1.0)));
        // Saving deltaS_R in nodeVariables
        nodeVariables[node][numVars] = deltaS_R[node];
    }

    // Function g(deltaS_R) on nodes --> Paparone - Tognaccini drag formulation
    vector <double> g_values(numNodes, 0.0);
    for (int node = 0; node < numNodes; node++) {
        g_values[node] = - fs1 * (deltaS_R[node]) - fs2 * pow((deltaS_R[node]), 2.0);
    }

    // Non dimensional enthalpy variation on nodes
    vector <double> NonDimDeltaH(numNodes, 0.0);
    for (int node = 0; node < numNodes; node++) {
        double T   = nodeVariables[node][T_id - 1]; // T on node
        double rho = nodeVariables[node][rho_id - 1]; // rho on node
        double u   = nodeVariables[node][momentumX_id - 1] / rho;
        double v   = nodeVariables[node][momentumY_id - 1] / rho;
        NonDimDeltaH[node] = ((R_air*gamma_air*T/(gamma_air - 1) + ((pow(u,2) + pow(v,2)) / 2)) - (gamma_air*R_air*T_inf / (gamma_air - 1) + ((pow(V_inf,2)) / 2))) / pow(V_inf,2);
        nodeVariables[node][numVars + 1] = NonDimDeltaH[node];
    }

    // Function f = V/V_inf = f(NonDimDeltaH) --> Destarac - Van der Vooren drag formulation
    vector <double> f_values(numNodes, 0.0);
    for (int node = 0; node < numNodes; node++) {
        f_values[node] = - sqrt(1 + 2.0 * (NonDimDeltaH[node]) - 2.0 / ((gamma_air - 1) * pow(M_inf,2.0)) * ((exp(deltaS_R[node] * (gamma_air - 1) / gamma_air)) - 1));
        if (std::isnan(f_values[node])) {
            f_values[node] = 0.0;
        }
    }

    /* Quantities to compute drag and lift in Lamb vector-based methods */

    // Cell centered values of vector r = X - X_pole --> Pole coordinates read from configuration file (.cas)
    double r_cc[2] = {0.0, 0.0};
    for (int cell = 0; cell < numCells; cell++) {
        double centroid_x = cellCenteredVariables[cell][CENTROID_X];
        double centroid_y = cellCenteredVariables[cell][CENTROID_Y];
        r_cc[0] = centroid_x - x_pole;
        r_cc[1] = centroid_y - y_pole;
        ccVariables[cell][R_XC] = r_cc[0];
        ccVariables[cell][R_YC] = r_cc[1];
    }

    // Nodal values of vector r = X - X_pole --> Pole coordinates read from configuration file (.cas)
    vector<double> r_x(numNodes, 0.0);
    vector<double> r_y(numNodes, 0.0);
    for (int node = 0; node < numNodes; node++) {
        double x = nodeVariables[node][0];
        double y = nodeVariables[node][1];
        r_x[node] = x - x_pole;
        r_y[node] = y - y_pole;
        nodeVariables[node][numVars + 4] = r_x[node];
        nodeVariables[node][numVars + 5] = r_y[node];
    }

    // Compute cell centered density, non dimensional entropy, total enthalpy, kinetic energy, and turbulent kinetic energy gradients (GREEN-GAUSS node based method)
    for(int cell = 0; cell < numCells; cell++) {

        double A_cell = cellCenteredVariables[cell][CELL_AREA]; // cell area
        double grad_rho[2] = {0.0, 0.0}; // Drho/Dx, Drho/Dy
        double grad_S[2] = {0.0, 0.0}; // dS/dx, dS/dy
        double grad_H[2] = {0.0, 0.0}; // dH/dx, dH/dy
        double grad_K[2] = {0.0, 0.0}; // (u*du/dx + v*dv/dx, u*du/dy + v*dv/dy) or (d(V^2/2)/dx + d(V^2/2)/dy)
        double grad_Kt[2] = {0.0, 0.0}; // dKt/dx, dKt/dy

        // Compute gradient contributions at face centeres 
        for(int localFace = 0; localFace < 4; localFace++) {

            // Face nodes index
            int64_t node1 = connectivity[cell][localFace];
            int64_t node2 = connectivity[cell][(localFace + 1) % 4];

            if (node1 == node2) {
                continue;
            }

            // Nodal values
            double rho1  = nodeVariables[node1 - 1][rho_id - 1];
            double rho2  = nodeVariables[node2 - 1][rho_id - 1];
            double S1    = nodeVariables[node1 - 1][numVars];
            double S2    = nodeVariables[node2 - 1][numVars];
            double H1    = nodeVariables[node1 - 1][numVars + 1];
            double H2    = nodeVariables[node2 - 1][numVars + 1];
            double Kt1   = nodeVariables[node1 - 1][k_id - 1];
            double Kt2   = nodeVariables[node2 - 1][k_id - 1];
            double u1    = nodeVariables[node1 - 1][momentumX_id - 1] / rho1;
            double u2    = nodeVariables[node2 - 1][momentumX_id - 1] / rho2;
            double v1    = nodeVariables[node1 - 1][momentumY_id - 1] / rho1;
            double v2    = nodeVariables[node2 - 1][momentumY_id - 1] / rho2;
            double K1    = 0.5 * (u1*u1 + v1*v1);
            double K2    = 0.5 * (u2*u2 + v2*v2);

            // Face centered interpolation
            double rho_face = 0.5 * (rho1 + rho2);
            double S_face = 0.5 * (S1 + S2);
            double H_face = 0.5 * (H1 + H2);
            double Kt_face = 0.5 * (Kt1 + Kt2);
            double u_face = 0.5 * (u1 + u2);
            double v_face = 0.5 * (v1 + v2);
            double K_face = 0.5 * (K1 + K2);

            // Normals
            double n_x = normals[cell][localFace * 2];
            double n_y = normals[cell][localFace * 2 + 1];

            // Sum on faces
            grad_rho[0] += rho_face * n_x;
            grad_rho[1] += rho_face * n_y;
            grad_S[0] += S_face * n_x;
            grad_S[1] += S_face * n_y;
            grad_H[0] += H_face * n_x;
            grad_H[1] += H_face * n_y;
            grad_Kt[0] += Kt_face * n_x;
            grad_Kt[1] += Kt_face * n_y;
            grad_K[0] += K_face * n_x;
            grad_K[1] += K_face * n_y;
        }

        // Cell centered gradient components
        grad_rho[0] /= A_cell;
        grad_rho[1] /= A_cell;
        grad_S[0] /= A_cell;
        grad_S[1] /= A_cell;
        grad_H[0] /= A_cell;
        grad_H[1] /= A_cell;
        grad_Kt[0] /= A_cell;
        grad_Kt[1] /= A_cell;
        grad_K[0] /= A_cell;
        grad_K[1] /= A_cell;

        // Multiply for R because S is the non-dimensional enthropy variation
        grad_S[0] *= R_air;
        grad_S[1] *= R_air;
        // Multiply for V_inf^2 because H is the non-dimensional enthalpy variation
        grad_H[0] *= pow(V_inf,2);
        grad_H[1] *= pow(V_inf,2);

        // // Cell centered terms to compute dK/dx and dK/dy --> IN THIS WAY A NUMERICAL ERROR IS INTRODUCED! => REPLACED
        // double u_cc  = ccVariables[cell][U_CC];
        // double v_cc  = ccVariables[cell][V_CC];
        // double V_magn = sqrt (u_cc*u_cc + v_cc*v_cc);
        // double du_dx = cellCenteredVariables[cell][DU_DX];
        // double du_dy = cellCenteredVariables[cell][DU_DY];
        // double dv_dx = cellCenteredVariables[cell][DV_DX];
        // double dv_dy = cellCenteredVariables[cell][DV_DY];

        // // Compute grad_K 
        // grad_K[0] = u_cc * du_dx + v_cc * dv_dx;
        // grad_K[1] = u_cc * du_dy + v_cc * dv_dy;

        // Saving gradient components 
        cellCenteredVariables[cell][DRHO_DX] = grad_rho[0];  // drho/dx
        cellCenteredVariables[cell][DRHO_DY] = grad_rho[1];  // drho/dy 
        cellCenteredVariables[cell][DS_DX]   = grad_S[0];  // dS/dx
        cellCenteredVariables[cell][DS_DY]   = grad_S[1];  // dS/dy
        cellCenteredVariables[cell][DH_DX]   = grad_H[0];  // dH/dx
        cellCenteredVariables[cell][DH_DY]   = grad_H[1];  // dH/dy
        cellCenteredVariables[cell][DK_DX]   = grad_K[0];  // dK/dx
        cellCenteredVariables[cell][DK_DY]   = grad_K[1];  // dK/dy
        cellCenteredVariables[cell][DKT_DX]  = grad_Kt[0];  // dKt/dx
        cellCenteredVariables[cell][DKT_DY]  = grad_Kt[1];  // dKt/dy
    }

    // WARNING: IN TAU_V MUST BE CONSIDERED ONLY mu_d AND NOT mu_t!
    // Compute again viscous stress tensor and its divergence considering just mu_d for far field calculation
    for (int cell = 0; cell < numCells; cell++) {
        double numFaces = 0.0;

        // Velocity gradient components
        double du_dx = cellCenteredVariables[cell][DU_DX];
        double du_dy = cellCenteredVariables[cell][DU_DY];
        double dv_dx = cellCenteredVariables[cell][DV_DX];
        double dv_dy = cellCenteredVariables[cell][DV_DY];

        // Cell centered dynamic viscosity
        double mu_d_cc = ccVariables[cell][MU_D_CC];

        // Cell centered velocity field divergence 
        double divV = du_dx + dv_dy; // du/dx + dv/dy

        // Viscous stress tensor components
        double tau11  = 2 * mu_d_cc * (du_dx - (divV / 3)); // 2 * mu * (du/dx - divV/3)
        double tau12  = mu_d_cc * (dv_dx + du_dy);      // 2 * mu * ((dv/dx + du/dy)/2)
        //  double tau21 = tau12; (simmetry)
        double tau22  = 2 * mu_d_cc * (dv_dy - (divV / 3));  // 2 * mu * (dv/dy - divV/3) 

        // Saving tau components in cellCenteredVariables matrix
        cellCenteredVariables[cell][TAU11_FF] = tau11;
        cellCenteredVariables[cell][TAU12_FF] = tau12; 
        cellCenteredVariables[cell][TAU21_FF] = tau12;
        cellCenteredVariables[cell][TAU22_FF] = tau22;
    }

    // Compute components of tau gradient with Green-Gauss CELL BASED (tau_f = a*tau_p + (1-a)*tau_n) considering only mu_d
    for (int cell = 0; cell < numCells; cell++) {

        double dTaudx[2] = {0.0, 0.0};
        double dTaudy[2] = {0.0, 0.0};
        double cellArea = cellCenteredVariables[cell][CELL_AREA];

        for (int localFace = 0; localFace < 4; localFace++) {
            int64_t node1 = connectivity[cell][localFace];
            int64_t node2 = connectivity[cell][(localFace + 1) % 4];

            if (node1 == node2) continue;

            double n_x = normals[cell][localFace * 2];
            double n_y = normals[cell][localFace * 2 + 1];

            // Face centroid
            double x_f = 0.5 * (nodeVariables[node1 - 1][0] + nodeVariables[node2 - 1][0]);
            double y_f = 0.5 * (nodeVariables[node1 - 1][1] + nodeVariables[node2 - 1][1]);

            // Looking for adjacent cell
            pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
            int adjacentCell = -1; // none
            auto it = faceToCells.find(face); // reference to faceToCells
            if (it != faceToCells.end() && it->second.size() == 2) { // check if 'face' is shared exactly by two cells
                adjacentCell = (it->second[0] == cell) ? it->second[1] : it->second[0]; // define adjacent cell (different from current cell)
                //cout << "FOUND N1: " << face.first << ", N2: " << face.second << ", WORKING CELL: " << cell << ", ADJACENT CELL: " << adjacentCell << "\n";
            }

            // Face centered tau interpolation
            double tau11_face = 0.0, tau12_face = 0.0, tau22_face = 0.0;

            if(adjacentCell != -1) { // if adjacent cell exists
                double tau11_n = cellCenteredVariables[cell][TAU11_FF];
                double tau11_p = cellCenteredVariables[adjacentCell][TAU11_FF];
                double tau12_n = cellCenteredVariables[cell][TAU12_FF];
                double tau12_p = cellCenteredVariables[adjacentCell][TAU12_FF];
                double tau22_n = cellCenteredVariables[cell][TAU22_FF];
                double tau22_p = cellCenteredVariables[adjacentCell][TAU22_FF];

                // Compute a 
                double a = 0.0;
                double dist_n = sqrt(pow(x_f - cellCenteredVariables[cell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[cell][CENTROID_Y], 2));
                double dist_p = sqrt(pow(x_f - cellCenteredVariables[adjacentCell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[adjacentCell][CENTROID_Y], 2));
                a = dist_n / (dist_p + dist_n);
                
                tau11_face = a * tau11_p + (1.0 - a) * tau11_n;
                tau12_face = a * tau12_p + (1.0 - a) * tau12_n;
                tau22_face = a * tau22_p + (1.0 - a) * tau22_n;
            } else {
                tau11_face = cellCenteredVariables[cell][TAU11_FF];
                tau12_face = cellCenteredVariables[cell][TAU12_FF];
                tau22_face = cellCenteredVariables[cell][TAU22_FF];
            }

            // Tau gradient components
            dTaudx[0] += tau11_face * n_x; // dtau11/dx
            dTaudy[1] += tau12_face * n_y; // dtau12/dy
            dTaudx[1] += tau12_face * n_x; // dtau12/dx
            dTaudy[0] += tau22_face * n_y; // dtau22/dy
        }

        dTaudx[0] /= cellArea; // dtau11/dx
        dTaudy[1] /= cellArea; // dtau12/dy
        dTaudx[1] /= cellArea; // dtau12/dx
        dTaudy[0] /= cellArea; // dtau22/dy
        
        ccVariables[cell][DTAU11DX_FF] = dTaudx[0];
        ccVariables[cell][DTAU12DY_FF] = dTaudy[1];
        ccVariables[cell][DTAU12DX_FF] = dTaudx[1];
        ccVariables[cell][DTAU22DY_FF] = dTaudy[0];
    }

    // Compute Lamb vector and integrand quantities included in Fs volume integral
    // 1) Compute cell-centered quantities: rho(r dot l) and rho*r*l tensor, where r is the position vector, l is the lamb vector to be computed as omega x V in std approach and T*gradS - gradH + div(tau_v)/rho in thd approach
    for (int cell = 0; cell < numCells; cell++) {
        double l_std_x = 0.0;
        double l_std_y = 0.0;
        double l_thd_x = 0.0;
        double l_thd_y = 0.0;
        double r_dot_l_std = 0.0, r_dot_l_thd = 0.0;
        double r_rhol_std[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
        double r_rhol_thd[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
        double centroid_x = cellCenteredVariables[cell][CENTROID_X];
        double centroid_y = cellCenteredVariables[cell][CENTROID_Y];
        double r_x = ccVariables[cell][R_XC];
        double r_y = ccVariables[cell][R_YC];
        double rho_cc = ccVariables[cell][RHO_CC];

        // Compute Lamb vector with standard approach (omega x V)
        if(vrt_std_on == 1) { 
            double omega = cellCenteredVariables[cell][DV_DX] - cellCenteredVariables[cell][DU_DY];
            l_std_x = - omega * ccVariables[cell][V_CC];
            l_std_y =   omega * ccVariables[cell][U_CC];
        }

        // Compute Lamb vector wit thermodynamic approach (T*gradS - gradH + div(tau_v)/rho)
        double T_cc = ccVariables[cell][T_CC];
        double dSdx = cellCenteredVariables[cell][DS_DX];
        double dSdy = cellCenteredVariables[cell][DS_DY];
        double dHdx = cellCenteredVariables[cell][DH_DX];
        double dHdy = cellCenteredVariables[cell][DH_DY];
        double divTau_x = ccVariables[cell][DTAU11DX_FF] + ccVariables[cell][DTAU12DY_FF];
        double divTau_y = ccVariables[cell][DTAU12DX_FF] + ccVariables[cell][DTAU22DY_FF];
        l_thd_x = (T_cc * dSdx - dHdx) + (divTau_x/rho_cc);
        l_thd_y = (T_cc * dSdy - dHdy) + (divTau_y/rho_cc);

        // Compute rho(r dot l) for both approaches
        r_dot_l_std = (r_x * l_std_x + r_y * l_std_y) * rho_cc; //std
        r_dot_l_thd = (r_x * l_thd_x + r_y * l_thd_y) * rho_cc; //thd

        // Compute rho*r*l tensor for both approaches
        r_rhol_std[0][0] = r_x * rho_cc * l_std_x;
        r_rhol_std[0][1] = r_x * rho_cc * l_std_y;
        r_rhol_std[1][0] = r_y * rho_cc * l_std_x;
        r_rhol_std[1][1] = r_y * rho_cc * l_std_y;
        r_rhol_thd[0][0] = r_x * rho_cc * l_thd_x;
        r_rhol_thd[0][1] = r_x * rho_cc * l_thd_y;
        r_rhol_thd[1][0] = r_y * rho_cc * l_thd_x;
        r_rhol_thd[1][1] = r_y * rho_cc * l_thd_y;

        // Store computed values in auxiliary structures
        cellCenteredVariables[cell][LAMB_STD_X] = l_std_x;
        cellCenteredVariables[cell][LAMB_STD_Y] = l_std_y;
        cellCenteredVariables[cell][LAMB_THD_X] = l_thd_x;
        cellCenteredVariables[cell][LAMB_THD_Y] = l_thd_y;
        ccVariables[cell][R_DOT_L_std] = r_dot_l_std;
        ccVariables[cell][R_DOT_L_thd] = r_dot_l_thd; 
        ccVariables[cell][R_RHOL_std_XX] = r_rhol_std[0][0];
        ccVariables[cell][R_RHOL_std_XY] = r_rhol_std[0][1];
        ccVariables[cell][R_RHOL_std_YX] = r_rhol_std[1][0];
        ccVariables[cell][R_RHOL_std_YY] = r_rhol_std[1][1];
        ccVariables[cell][R_RHOL_thd_XX] = r_rhol_thd[0][0];
        ccVariables[cell][R_RHOL_thd_XY] = r_rhol_thd[0][1];
        ccVariables[cell][R_RHOL_thd_YX] = r_rhol_thd[1][0];
        ccVariables[cell][R_RHOL_thd_YY] = r_rhol_thd[1][1];
    }

    // ++++++++++++++++++++++++++ PART 2:VARIABLE QUANTITIES WITH CONTROL VOLUME => LOOP ON X_wake ++++++++++++++++++++++++++
    // +++++++++++++++++++++ EXTERNAL LOOP ON X_W TO DO CALCULATION SELECTING DIFFERENT CONTROL VOLUMES +++++++++++++++++++++

    ofstream outputFile("PPUN2D_OUT.txt");
    if (vrt_std_on == 1) {
        outputFile << "Xmax\t\tCd_PT\t\tCd_PT_wv\t\tCd_PT_vs\t\tCd_PT_sp\t\tCd_DV\t\tCd_DV_wv\t\tCd_DV_vs\t\tCd_DV_sp\t\tCd_vf_thd\t\tCl_vf_thd\t\tCd_vf_std\t\tCl_vf_std\t\tCd_compr\t\tCl_compr\t\tCd_fs_thd\t\tCl_fs_thd\t\tCd_fs_std\t\tCl_fs_std\t\tCd_fmu_thd\t\tCl_fmu_thd\t\tCd_fmu_std\t\tCl_fmu_std\t\tCd_vrt_thd\t\tCl_vrt_thd\t\tCd_vrt_std\t\tCl_vrt_std\t\tCd_fs_vol_thd\t\tCd_fs_vol_std\t\tCd_p_thd\t\tCd_p_std\t\tCd_i_thd\t\tCd_i_std\t\tCd_wv_thd\t\tCd_vs_thd\t\tCd_sp_thd\t\tCd_wv_std\t\tCd_vs_std\t\tCd_sp_std\t\tCd_nf\t\tCl_nf\n";
    } else {
        outputFile << "Xmax\t\tCd_PT\t\tCd_PT_wv\t\tCd_PT_vs\t\tCd_PT_sp\t\tCd_DV\t\tCd_DV_wv\t\tCd_DV_vs\t\tCd_DV_sp\t\tCd_vf_thd\t\tCl_vf_thd\t\tCd_compr\t\tCl_compr\t\tCd_fs_thd\t\tCl_fs_thd\t\tCd_fmu\t\tCl_fmu\t\tCd_vrt_thd\t\tCl_vrt_thd\t\tCd_fs_vol_thd\t\tCd_p_thd\t\tCd_i_thd\t\tCd_wv_thd\t\tCd_vs_thd\t\tCd_sp_thd\t\tCd_nf\t\tCl_nf\n";
    }

    double elapsed_time_single = 0.0; // execution time of one iteration + previous calculation
    
    double x_max_i = 1.0, x_max_step = 1;
    for (double Xmax = x_max_i; Xmax <= x_max; Xmax += x_max_step) {

        // Variables Initialization
        // Far field (THD)
        double I_D_PT_tot = 0.0; // integrand Paparone - Tognaccini
        double I_D_PT_wave = 0.0;
        double I_D_PT_viscous = 0.0;
        double I_D_PT_spurious = 0.0;
        double D_PT_tot = 0.0, D_PT_wave = 0.0, D_PT_viscous = 0.0, D_PT_spurious = 0.0; // entropy drag Paparone - Tognaccini
        double I_D_DV_tot = 0.0; // integrand Destarac - Van der Vooren
        double I_D_DV_wave = 0.0;
        double I_D_DV_viscous = 0.0;
        double I_D_DV_spurious = 0.0;
        double D_DV_tot = 0.0, D_DV_wave = 0.0, D_DV_viscous = 0.0, D_DV_spurious = 0.0; // entropy drag Paparone - Tognaccini
        // Far field (Lamb Vector)
        double F_l_x = 0.0;
        double F_l_y = 0.0;
        double Fl_std_x = 0.0;
        double Fl_std_y = 0.0;
        double F_s_x = 0.0; // surface integral (thd approach)
        double F_s_y = 0.0; // surface integral (thd approach)
        double Fs_std_x = 0.0; // surface integral
        double Fs_std_y = 0.0; // surface integral
        double F_mrho_x = 0.0;
        double F_mrho_y = 0.0;
        double F_mud_x_thd = 0.0;
        double F_mud_y_thd = 0.0;
        double F_mu_x_std = 0.0;
        double F_mu_y_std = 0.0;
        double F_mut_x = 0.0;
        double F_mut_y = 0.0;
        double F_mu_x_thd = 0.0;
        double F_mu_y_thd = 0.0;
        double F_mut2_x = 0.0, F_mut2_y = 0.0;
        double F_mut1_x = 0.0, F_mut1_y = 0.0;
        // NodeID flag initialization
        for (int node = 0; node < numNodes; node++) {
            nodeVariables[node][numVars + 2] = 0.0;
            nodeVariables[node][numVars + 3] = 0.0;
        }

        /* Definition of the internal (body) and external (box) boundary of the control volume */
        // Creation of map 'Face counter'
        map < pair < int64_t, int64_t>, int> faceCount;
        faceCount.clear();
        // Selecting cells and faces inside the control volume
        vector <int64_t> cellsInControlVolume;
        cellsInControlVolume.clear();
        for (int64_t cell = 0; cell < numCells; cell++) {
            
            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];

            if(centroid_x >= x_min && centroid_x <= Xmax  && centroid_y >= y_min && centroid_y <= y_max) {
                cellsInControlVolume.push_back(cell);

            int numFaces = 4;
                // Loop on faces to update faceCount
                for(int localFace = 0; localFace < numFaces; localFace++) {
                    int64_t node1 = connectivity[cell][localFace]; // first node of localFace
                    int64_t node2 = connectivity[cell][(localFace + 1) % numFaces]; // second node of localFace
                    if(node1 == node2) {
                        continue;
                    } 
                    // Node sorting to avoid duplicates (FUNDAMENTAL)
                    pair <int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                    // Update face counter 
                    faceCount[face]++; // if 'face' (key) doesn't exist, then its value is equal to 1. If face exists then its value is updated (1 + 1 + ...)
                }
            }
        }

        // Selecting nodes on boundaries
        for (const pair<pair<int64_t, int64_t>, int>& entry : faceCount) {
            const pair<int64_t, int64_t>& face = entry.first; // 'face' (key)
            int count = entry.second; // 'count' (value)
            // Check if face is on boundaries
            if (count == 1) { // face not shared with oter cells => faces on boundary
                int64_t node1 = face.first;
                int64_t node2 = face.second;

                // Nodes coordinates
                double x1 = nodeVariables[node1 - 1][0];
                double y1 = nodeVariables[node1 - 1][1];
                double x2 = nodeVariables[node2 - 1][0];
                double y2 = nodeVariables[node2 - 1][1];
    
                // Normals
                double n_x = y2 - y1;
                double n_y = -(x2 - x1);

                // Check on Vn to select faces on external or internal boundary
                double V_x1 = nodeVariables[node1 - 1][momentumX_id - 1] / nodeVariables[node1 - 1][rho_id - 1];
                double V_y1 = nodeVariables[node1 - 1][momentumY_id - 1] / nodeVariables[node1 - 1][rho_id - 1];
                double V_x2 = nodeVariables[node2 - 1][momentumX_id - 1] / nodeVariables[node2 - 1][rho_id - 1];
                double V_y2 = nodeVariables[node2 - 1][momentumY_id - 1] / nodeVariables[node2 - 1][rho_id - 1];

                double V_n1 = V_x1 * n_x + V_y1 * n_y;
                double V_n2 = V_x2 * n_x + V_y2 * n_y;

                if (fabs(V_n1 < 1e-6 && fabs(V_n2) < 1e-6)) {
                    nodeVariables[node1 - 1][numVars + 3] = 1; // nodeID_body
                    nodeVariables[node2 - 1][numVars + 3] = 1; // nodeID_body
                } else { 
                    nodeVariables[node1 - 1][numVars + 2] = 1; // nodeID_ext
                    nodeVariables[node2 - 1][numVars + 2] = 1; // nodeID_ext
                }
            }
        }

        /* ----------- NEAR-FIELD Aerodynamic force calculation --> D_near = Integral[((p-p_inf)n - tau*n)dS] ---------- */
        double Fx_nf = 0.0;
        double Fy_nf = 0.0;
        double Fx_p = 0.0, Fx_f, Fy_p = 0.0, Fy_f = 0.0;
        if (Xmax == x_max_i) { // Near-field drag is not dipendent on Xmax 
            for (int64_t cell = 0; cell < numCells; cell++) { 
                for (int localFace = 0; localFace < 4; localFace++) {
                    // Face nodes index
                    int64_t node1 = connectivity[cell][localFace];
                    int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                    if (node1 == node2) {
                        continue;
                    }

                    // Check if nodes are on body surface
                    if (nodeVariables[node1 - 1][numVars + 3] == 1 && nodeVariables[node2 - 1][numVars + 3] == 1) {
                        
                    // Face centered pressure interpolation
                    double p1 = nodeVariables[node1 - 1][p_id - 1];
                    double p2 = nodeVariables[node2 - 1][p_id - 1];
                    double p_face = 0.5 * (p1 + p2);
                

                    // Normals 
                    double n_x = normals[cell][localFace * 2];
                    double n_y = normals[cell][localFace * 2 + 1];

                    // Pressure contribution
                    double pn_x = p_face * n_x;
                    double pn_y = p_face * n_y;
                    Fx_p += pn_x; // Cd_p 
                    Fy_p += pn_y; // Cl_p 

                    // Viscous contribution
                    double tau11 = cellCenteredVariables[cell][TAU11_NF];
                    double tau12 = cellCenteredVariables[cell][TAU12_NF];
                    double tau22 = cellCenteredVariables[cell][TAU22_NF];
                    double tn_x = tau11 * n_x + tau12 * n_y;
                    double tn_y = tau12 * n_x + tau22 * n_y;
                    Fx_f += - tn_x; // Cd_f
                    Fy_f += - tn_y; // Cl_f 

                    // Compute total force contribution for localCell and sum (body normals are oppsite then external boundary normals)
                    Fx_nf += (pn_x - tn_x);
                    Fy_nf += (pn_y - tn_y); 
                    }
                }
            } 

            // PRINT NEAR-FIELD CALCULATION 
            std::cout << "-----------NEAR-FIELD METHOD RESULTS-----------\n";
            std::cout << "\n";
            std::cout << "AERODYNAMIC FORCE CONTRIBUTIONS:\n"
                << "Cd_p = " << Fx_p / (q_inf * chord) << "; Cd_f = " << Fx_f / (q_inf * chord) << "\n"
                << "Cl_p = " << Fy_p / (q_inf * chord) << ";  Cl_f = " << Fy_f / (q_inf * chord) << endl;
            std::cout << "\n";
            std::cout << "TOTAL AERODYNAMIC FORCE:\n"
                << "Drag (Fx) = " << Fx_nf << " N \n"
                << "Lift (Fy) = " << Fy_nf << " N" << endl;
            std::cout << "\n";
            std::cout << "AERODYNAMIC FORCE COEFFICIENTS:\n"
                << "Drag coefficient: Cd,nf = " << Fx_nf / (q_inf * chord) << "\n"
                << "Lift coefficient: Cl,nf = " << Fy_nf / (q_inf * chord) << "\n" << endl;
            if(split == 1) {
                std::cout << "WARNING: IF SPLITTING IS ON Cd,nf IS UNDERESTIMATED!" << endl;
            }

            // Compute Cd anc Cl only for the 1st iteration, then consider always the same constant value
            Cd_nf_stored = Fx_nf / (q_inf * chord);
            Cl_nf_stored = Fy_nf / (q_inf * chord);
        }

        // Put Cd_nf and Cl_nf = stored values for X_w > 1 
        double Cd_nf = Cd_nf_stored;
        double Cl_nf = Cl_nf_stored;

        /* Sensors to select viscous, shock and spurious regions */

        // Cell with max distance from body (0,0)
        double d_max = 0.0;
        int farthestCell = -1;
        for (int cell = 0; cell < numCells; cell++) {
            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];
            double d = sqrt(centroid_x*centroid_x + centroid_y*centroid_y);
            if(d > d_max) {
                d_max = d;
                farthestCell = cell;
            }
        }

        // F_bl_inf (Free stream value of F_bl)
        double F_bl_inf = (ccVariables[farthestCell][MU_D_CC] + ccVariables[farthestCell][MU_T_CC]) / ccVariables[farthestCell][MU_D_CC];
        //std::cout << "Free stream value of F_bl: " << F_bl_inf << endl;

        // Flag to identify shock, viscous and spurious regions
        std::vector<int> regionFlags(numCells, 0.0); // 0 = spurious, 1 = shock, 2 = viscous

        // Loop to compute sensors at CELL CENTERS and to define shock and viscous regions
        for (int cell = 0; cell < numCells; cell++) {

            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];
            double u_cc = 0.0, v_cc = 0.0, T_cc = 0.0, rho_cc = 0.0, mu_d_cc = 0.0, mu_t_cc = 0.0, k_cc = 0.0, dp_dx = 0.0, dp_dy = 0.0;
        
            if (centroid_x > x_min && centroid_x < Xmax && centroid_y > y_min && centroid_y < y_max) {

                    // Cell centered variables
                    u_cc = ccVariables[cell][U_CC];
                    v_cc = ccVariables[cell][V_CC];
                    T_cc = ccVariables[cell][T_CC];
                    rho_cc = ccVariables[cell][RHO_CC];
                    mu_d_cc = ccVariables[cell][MU_D_CC];
                    mu_t_cc = ccVariables[cell][MU_T_CC];
                    k_cc = ccVariables[cell][K_CC];
                    dp_dx = ccVariables[cell][DP_DX];
                    dp_dy = ccVariables[cell][DP_DY];

                    // Function F_shock  = (V dot gradP) / (a*|gradP|) for shock region
                    double F_shock = (u_cc * dp_dx + v_cc * dp_dy) / (sqrt(gamma_air * R_air * T_cc) * sqrt(dp_dx*dp_dx + dp_dy*dp_dy));

                    // Function F_bl = mu_d + mu_t / mu_d for viscous region
                    double F_bl = (mu_d_cc + mu_t_cc) / mu_d_cc;

                    // Omega (turbulence dissipation frequency) BL sensor 
                    double OmegaSensor = rho_cc * k_cc / (std::max(mu_t_cc, 1e-15));

                    // Region selection
                    if (F_shock > K_w) {
                        regionFlags[cell] = 1; // shock wave region
                    }
                    else if (OmegaSensor > omega_tsh || F_bl > K_bl * F_bl_inf) {
                        regionFlags[cell] = 2; // boundary layer and wake region
                    }
                    else {
                        regionFlags[cell] = 0; // spurious region
                    }

                    cellCenteredVariables[cell][CELL_TYPE] = regionFlags[cell];
            }
        }

        // Expanding the shock wave region. In this way drag oscillations near the shock wave are taking into account.
        
        // Step 1: Identify cells adjacent currently included in the shock wave region
        unordered_set<int64_t> shockCells;
        for (int cell = 0; cell < numCells; cell++) {
            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];
            if (centroid_x > x_min && centroid_x < Xmax && centroid_y > y_min && centroid_y < y_max) {
                if (regionFlags[cell] == 1) {  // If cell is currently identified as a shock one
                    shockCells.insert(cell); // Store for further propagation     
                }
            }
        }
            
        // Step 2: Expand shock wave region outward for the first X = N_shock_margin - 1 layers
        unordered_set<int64_t> newShockLayer = shockCells;

        for (int i = 0; i < N_shock_margin; i++) {
            unordered_set<int64_t> nextLayer;
            
            for (int64_t cell : newShockLayer) {
                for (int localFace = 0; localFace < 4; localFace++) {
                    int64_t node1 = connectivity[cell][localFace];
                    int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                    // Get adjacent cells
                    pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                    auto it = faceToCells.find(face);
                    if (it != faceToCells.end()) {
                        for (int adjCell : it->second) {
                            if (adjCell != cell && regionFlags[adjCell] == 0) { // Expand only to unmarked cells
                                regionFlags[adjCell] = 1; // Set as shock
                                cellCenteredVariables[adjCell][CELL_TYPE] = 1;
                                nextLayer.insert(adjCell);
                            }
                        }
                    }
                }
            }
            newShockLayer = nextLayer; // Move to the next layer
        } 

        // Expand the shock wave region in streamwise direction to consider the shock wake, if any
        unordered_set<int64_t> shockBoundaryCells;
        double max_x_shock = 1e-9;

        for (int64_t cell = 0; cell < numCells; cell++) { // find the maximum x among shock cells
            if (regionFlags[cell] == 1) { // if cell is a shock one
                double x = cellCenteredVariables[cell][CENTROID_X];
                if (x > max_x_shock) {
                    max_x_shock = x;
                }
            }
        }

        double min_y_shock = 1e9;
        double max_y_shock = 1e-9;

        for (int64_t cell = 0; cell < numCells; cell++) { // find the maximum x among shock cells
            if (regionFlags[cell] == 1) { // if cell is a shock one
                double y = cellCenteredVariables[cell][CENTROID_Y];
                if (y < min_y_shock) min_y_shock = y;
                if (y > max_y_shock) max_y_shock = y;
            }
        }

        // Collect shock cells near max_x_shock
        double x_tol = 0.05;
        for (int64_t cell = 0; cell < numCells; cell++) {
            if (regionFlags[cell] == 1) {
                double x = cellCenteredVariables[cell][CENTROID_X];
                if (x > max_x_shock - x_tol) {
                    shockBoundaryCells.insert(cell);
                }
            }
        }

        // Extend shock region downstream
        unordered_set<int64_t> currentShockLayer = shockBoundaryCells;

        for (int64_t layer = 0; layer < N_shock_streamwise_margin; layer++) {
            unordered_set<int64_t> nextLayer;

            for (int64_t cell : currentShockLayer) {
                double x = cellCenteredVariables[cell][CENTROID_X];

                for (int localFace = 0; localFace < 4; localFace++) {
                    int64_t node1 = connectivity[cell][localFace];
                    int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                    // Get adjacent cells
                    pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                    auto it = faceToCells.find(face);
                    if (it != faceToCells.end()) {
                        for (int adjCell : it->second) {
                            double x_adj = cellCenteredVariables[adjCell][CENTROID_X];
                            if (adjCell != cell && (regionFlags[adjCell] == 0) && x_adj > x) { // Expand downstream only
                                double y_adj = cellCenteredVariables[adjCell][CENTROID_Y];
                                if (y_adj >= min_y_shock && y_adj <= max_y_shock && x_adj <= x_max) {
                                    regionFlags[adjCell] = 1; // Remap to shock
                                    cellCenteredVariables[adjCell][CELL_TYPE] = 1;
                                    nextLayer.insert(adjCell);
                                }
                            }
                        }
                    }
                }
            }
            currentShockLayer = nextLayer;
        }

        // Expanding the viscous region
        // The actual region selection does not set the cells near the body as viscous, but as spurious. This is incorrect.
        // It is necessary to set the cells near the body (boundary layer) as viscous, including the first X layers outward.

        // Step 1: Identify cells adjacent to the body and set them as viscous
        unordered_set<int64_t> viscousCells;
        for (int cell = 0; cell < numCells; cell++) {
            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];

            if (centroid_x > x_min && centroid_x < Xmax && centroid_y > y_min && centroid_y < y_max) {
                if (regionFlags[cell] == 0) {  // If cell is currently spurious
                    for (int localFace = 0; localFace < 4; localFace++) {
                        int64_t node1 = connectivity[cell][localFace];
                        int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                        if (node1 == node2) continue; 

                        if (nodeVariables[node1 - 1][numVars + 3] == 1 && nodeVariables[node2 - 1][numVars + 3] == 1) { // If NodeID_body = true
                            regionFlags[cell] = 2; // Set as viscous region
                            cellCenteredVariables[cell][CELL_TYPE] = regionFlags[cell];
                            viscousCells.insert(cell); // Store for further propagation
                            break; // No need to check other faces
                        }
                    }
                }
            }
        }

        // Step 2: Expand viscous region outward for the first X = N_wall_margin - 1 layers
        unordered_set<int64_t> newViscousLayer = viscousCells;

        for (int i = 0; i < N_wall_margin; i++) {
            unordered_set<int64_t> nextLayer;
            
            for (int64_t cell : newViscousLayer) {
                for (int localFace = 0; localFace < 4; localFace++) {
                    int64_t node1 = connectivity[cell][localFace];
                    int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                    // Get adjacent cells
                    pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                    auto it = faceToCells.find(face);
                    if (it != faceToCells.end()) {
                        for (int adjCell : it->second) {
                            if (adjCell != cell && regionFlags[adjCell] == 0) { // Expand only to unmarked cells
                                regionFlags[adjCell] = 2; // Set as viscous
                                cellCenteredVariables[adjCell][CELL_TYPE] = regionFlags[adjCell];
                                nextLayer.insert(adjCell);
                            }
                        }
                    }
                }
            }
            newViscousLayer = nextLayer; // Move to the next layer
        }

        // Expanding the viscous region outward for a better selection of the viscous region.
        // The BL region on the body is expanded with a fixed number of layers, while the wake region growth linearly in x-direction.

        // Step 0: find TE
        double TE = 0.0;
        for(int64_t cell = 0; cell < numCells; cell++){
            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            bool isOnBody = false;
            for (int localFace = 0; localFace < 4; localFace++) {
                int64_t node1 = connectivity[cell][localFace];
                int64_t node2 = connectivity[cell][(localFace + 1) % 4];
                if(nodeVariables[node1 - 1][numVars + 3] == 1 && nodeVariables[node2 - 1][numVars + 3] == 1){
                    isOnBody = true;
                    break;
                }
            }

            if(isOnBody && centroid_x > TE){
                TE = centroid_x;
            }
        } 
        //std::cout<<"TE: " << TE << "\n";

        // Step 1: identify boundary cells between viscous and spurious region
        unordered_set<int64_t> viscousBoundaryCells, wakeBoundaryCells;
        for(int64_t cell = 0; cell < numCells; cell++){
            if(cellCenteredVariables[cell][CELL_TYPE] == 2){ // if cell is viscous
                bool isBoundary = false;
                for (int localFace = 0; localFace < 4; localFace++) {
                    int64_t node1 = connectivity[cell][localFace];
                    int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                    // Get adjacent cells
                    pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                    auto it = faceToCells.find(face);
                    if (it != faceToCells.end()) {
                        for (int adjCell : it->second) {
                            if (cellCenteredVariables[adjCell][CELL_TYPE] == 0) { // If adjacent cell is spurious then cell is a viscous boundary cell
                                isBoundary = true;
                                break;
                            }
                        }
                    }
                }

                // Distinction between boundary cells around the body and in the wake
                double centroid_x = cellCenteredVariables[cell][CENTROID_X];
                if(isBoundary){
                    if(centroid_x < TE){
                        viscousBoundaryCells.insert(cell); // Cells around the body
                    } else {
                        wakeBoundaryCells.insert(cell); // Cells in the viscous wake
                    }
                }
            }
        }

        // Step 2: Expand the viscous region around the body with a fixed number of layers (N_Bl_margin)
        unordered_set<int64_t> newViscousBodyLayer = viscousBoundaryCells;
        for (int i = 0; i < N_Bl_margin; i++) {
            unordered_set<int64_t> nextLayer;
            for (int64_t cell : newViscousBodyLayer) {
                for (int localFace = 0; localFace < 4; localFace++) {
                    int64_t node1 = connectivity[cell][localFace];
                    int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                    // Get adjacent cells
                    pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                    auto it = faceToCells.find(face);
                    if (it != faceToCells.end()) {
                        for (int adjCell : it->second) {
                            double adj_x = cellCenteredVariables[adjCell][CENTROID_X];
                            if (adjCell != cell && regionFlags[adjCell] == 0) { // Expand only to unmarked cells
                                regionFlags[adjCell] = 2; // Set as viscous
                                cellCenteredVariables[adjCell][CELL_TYPE] = regionFlags[adjCell];
                                nextLayer.insert(adjCell);
                            }
                        }
                    }
                }
            }
            newViscousBodyLayer = nextLayer; // Move to the next layer
        }

        // Step 3: Re-identify boundary cells in the wake after body extension
        unordered_set<int64_t> extendedWakeBoundaryCells;
        for(int64_t cell = 0; cell < numCells; cell++){
            if(cellCenteredVariables[cell][CELL_TYPE] == 2){ // if cell is viscous
                double centroid_x = cellCenteredVariables[cell][CENTROID_X];
                if(centroid_x > TE && centroid_x < Xmax){ // only wake
                    for (int localFace = 0; localFace < 4; localFace++) {
                        int64_t node1 = connectivity[cell][localFace];
                        int64_t node2 = connectivity[cell][(localFace + 1) % 4];
    
                        // Get adjacent cells
                        pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                        auto it = faceToCells.find(face);
                        if (it != faceToCells.end()) {
                            for (int adjCell : it->second) {
                                if (adjCell != cell && cellCenteredVariables[adjCell][CELL_TYPE] == 0) { // If adjacent cell is spurious then cell is a viscous boundary cell
                                    extendedWakeBoundaryCells.insert(cell);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Step 3.1: Find the upper and lower bounds of the viscous region on the body for a smooth transition
        double max_y_bl = -1e9, min_y_bl = 1e9;
        bool found = false;
        for(int64_t cell = 0; cell < numCells; cell++){
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];
            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            if(centroid_x >= TE - 0.1 && centroid_x <= TE + 0.1 && regionFlags[cell] == 2){
                //std::cout << "Cell: " << cell << "--> x = " << centroid_x << ", y = " << centroid_y << "\n";
                max_y_bl = std::max(max_y_bl, centroid_y);
                min_y_bl = std::min(min_y_bl, centroid_y);
                found = true;
            }
        }

        if(!found){
            std::cerr << "WARNING: TE cells not found in last body layer - using userd defined y-range values\n";
            max_y_bl = max_y_bl_udf;
            min_y_bl = min_y_bl_udf;

        } /* else {
            std::cout<<"Found TE y-range from body extension: max_y_bl = " << max_y_bl << ", min_y_bl = " << min_y_bl << std::endl;
        } */

        //Step 4: Expand the viscous region in the wake with a progressive number of layers
        unordered_set<int64_t> currentWakeBoundary = extendedWakeBoundaryCells;
        for(int layer = 0; layer < N_wake_margin; layer++){
            unordered_set<int64_t> nextLayer;
            for(int64_t wakeCell : currentWakeBoundary){
                double centroid_x = cellCenteredVariables[wakeCell][CENTROID_X];
                double x_ratio = (centroid_x - TE) / (Xmax - TE);
                int N_layers = N_Bl_margin + (int)(x_ratio*wake_slope*(N_wake_margin - N_Bl_margin));
                N_layers = std::max(N_layers, N_Bl_margin);
                if(layer >= N_layers) continue;
                for (int localFace = 0; localFace < 4; localFace++) {
                    int64_t node1 = connectivity[wakeCell][localFace];
                    int64_t node2 = connectivity[wakeCell][(localFace + 1) % 4];

                    // Get adjacent cells
                    pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                    auto it = faceToCells.find(face);
                    if (it != faceToCells.end()) {
                        for (int adjCell : it->second) {
                            double x_adj = cellCenteredVariables[adjCell][CENTROID_X];
                            double y_adj = cellCenteredVariables[adjCell][CENTROID_Y];
                            if (adjCell != wakeCell && regionFlags[adjCell] == 0 && x_adj < Xmax && x_adj > x_min && y_adj < y_max && y_adj > y_min) { // Expand only to unmarked cells
                                if(abs(x_adj - TE) < 0.05 && (y_adj > max_y_bl || y_adj < min_y_bl)){ // Check on max and min BL region height for smooth transition 
                                    continue;
                                } 
                                regionFlags[adjCell] = 2; // Set as viscous
                                cellCenteredVariables[adjCell][CELL_TYPE] = regionFlags[adjCell];
                                nextLayer.insert(adjCell);
                            }
                        }
                    }
                }
            }
            currentWakeBoundary = nextLayer;
        }

        // WAVE-BL INTERACTION: Expand viscous region into the shock wave to model shock-boundary layer interaction

        // Find viscous cells at the interface with shock wave
        unordered_set<int64_t> viscousShockInterfaceCells;

        for (int64_t cell = 0; cell < numCells; cell++) {
            if (cellCenteredVariables[cell][CELL_TYPE] == 2) { // Cell is viscous
                for (int localFace = 0; localFace < 4; localFace++) {
                    int64_t node1 = connectivity[cell][localFace];
                    int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                    // Get adjacent cells
                    pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                    auto it = faceToCells.find(face);
                    if (it != faceToCells.end()) {
                        for (int adjCell : it->second) {
                            if (adjCell != cell && cellCenteredVariables[adjCell][CELL_TYPE] == 1) { // Adjacent cell is shock
                                viscousShockInterfaceCells.insert(cell);
                                break;
                            }
                        }
                    }
                }
            }
        }

        // Now expand viscous region into the shock region starting from interface cells
        unordered_set<int64_t> currentLayer = viscousShockInterfaceCells;

        for (int i = 0; i < N_wv_bl_interaction; i++) {
            unordered_set<int64_t> nextLayer;

            for (int64_t cell : currentLayer) {
                for (int localFace = 0; localFace < 4; localFace++) {
                    int64_t node1 = connectivity[cell][localFace];
                    int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                    // Get adjacent cells
                    pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                    auto it = faceToCells.find(face);
                    if (it != faceToCells.end()) {
                        for (int adjCell : it->second) {
                            if (adjCell != cell && regionFlags[adjCell] == 1) { // Expand only into shock cells
                                regionFlags[adjCell] = 2; // Remap to viscous
                                cellCenteredVariables[adjCell][CELL_TYPE] = 2;
                                nextLayer.insert(adjCell);
                            }
                        }
                    }
                }
            }

            currentLayer = nextLayer; // Proceed to next layer
        }

        // Loop to compute far-field drag: Paparone Tognaccini and Destarac-Van der Vooren methods
        for (int64_t cell = 0; cell < numCells; cell++) {

            double cell_flux_PT = 0.0; // cell flux PT (TOTAL)
            double cell_flux_DV = 0.0; // cell flux DV (TOTAL)

            double cell_flux_PT_shock = 0.0;
            double cell_flux_PT_viscous = 0.0;
            double cell_flux_PT_spurious = 0.0;

            double cell_flux_DV_shock = 0.0;
            double cell_flux_DV_viscous = 0.0;
            double cell_flux_DV_spurious = 0.0;

            // Check if 'cell' is inside the control volume
            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];
            if(centroid_x < x_min || centroid_x > Xmax || centroid_y < y_min || centroid_y > y_max) {
                continue;
            }
            else {
                for (int localFace = 0; localFace < 4; ++localFace) {

                    double face_flux_PT = 0.0; // face flux PT
                    double face_flux_DV = 0.0; // face flux DV
                    double rho_face = 0.0;
                    double momentum_x_face = 0.0;
                    double momentum_y_face = 0.0;
                    double g_face = 0.0;
                    double f_face = 0.0;

                    // Node index for each face
                    int64_t node1 = connectivity[cell][localFace];
                    int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                    if (node1 == node2) {
                        continue;
                    }

                    // Face centered interpolation of nodal values 
                    rho_face = (nodeVariables[node1 - 1][rho_id - 1] + nodeVariables[node2 - 1][rho_id - 1]) / 2.0;
                    g_face = (g_values[node1 - 1] + g_values[node2 - 1]) / 2.0; 
                    f_face = (f_values[node1 - 1] + f_values[node2 - 1]) / 2.0;
                    momentum_x_face = (nodeVariables[node1 - 1][momentumX_id - 1] + nodeVariables[node2 - 1][momentumX_id - 1]) / 2.0;
                    momentum_y_face = (nodeVariables[node1 - 1][momentumY_id - 1] + nodeVariables[node2 - 1][momentumY_id - 1]) / 2.0;

                    // Normals (n*dS)
                    double n_x = normals[cell][localFace * 2];
                    double n_y = normals[cell][localFace * 2 + 1];

                    // PT flux on face: (rho g V) · n*dS
                    face_flux_PT = g_face * (momentum_x_face * n_x + momentum_y_face * n_y);
                    // DV flux on face: (rho f V) · n*dS
                    face_flux_DV = f_face * (momentum_x_face * n_x + momentum_y_face * n_y);  

                    // Total
                    cell_flux_PT += face_flux_PT;
                    cell_flux_DV += face_flux_DV;

                    // DRAG BREAKDOWN: SHOCK, VISCOUS AND SPURIOUS CONTRIBUTIONS
                    if (cellCenteredVariables[cell][CELL_TYPE] == 1.0) { // Wave drag
                        cell_flux_PT_shock += face_flux_PT;
                        cell_flux_DV_shock += face_flux_DV;
                    } 
                    else if (cellCenteredVariables[cell][CELL_TYPE] == 2.0) { // Viscous drag
                        cell_flux_PT_viscous += face_flux_PT;
                        cell_flux_DV_viscous += face_flux_DV;
                    }
                    else if (cellCenteredVariables[cell][CELL_TYPE] == 0.0) { // Spurious drag
                        cell_flux_PT_spurious += face_flux_PT;
                        cell_flux_DV_spurious += face_flux_DV;
                    }
                }
            }

            // Saving flux value for each cell
            // Wave
            cellCenteredVariables[cell][CELL_FLUX_DRAG_PT_WV] = cell_flux_PT_shock;
            cellCenteredVariables[cell][CELL_FLUX_DRAG_DV_WV] = cell_flux_DV_shock;
            // Viscous
            cellCenteredVariables[cell][CELL_FLUX_DRAG_PT_VS] = cell_flux_PT_viscous;
            cellCenteredVariables[cell][CELL_FLUX_DRAG_DV_VS] = cell_flux_DV_viscous;
            // Spurious
            cellCenteredVariables[cell][CELL_FLUX_DRAG_PT_SP] = cell_flux_PT_spurious;
            cellCenteredVariables[cell][CELL_FLUX_DRAG_DV_SP] = cell_flux_DV_spurious;
            // Total value
            cellCenteredVariables[cell][CELL_FLUX_DRAG_PT] = cell_flux_PT; 
            cellCenteredVariables[cell][CELL_FLUX_DRAG_DV] = cell_flux_DV;

            // Cell centered value of div(tau_x) = dTau11/dx + dTau12/dy --> correction term for DV formula
            double divTau_x = ccVariables[cell][DTAU11DX_NF] + ccVariables[cell][DTAU12DY_NF];
            
            // Sum cell_flux (compute integral)
            if (cellCenteredVariables[cell][CELL_TYPE] == 1.0) { // Wave
                I_D_PT_wave += cell_flux_PT_shock; 
                I_D_DV_wave += cell_flux_DV_shock;
            } else if (cellCenteredVariables[cell][CELL_TYPE] == 2.0) { // Viscous
                I_D_PT_viscous += cell_flux_PT_viscous; 
                I_D_DV_viscous += cell_flux_DV_viscous;  
            } else if (cellCenteredVariables[cell][CELL_TYPE] == 0.0) { // Spurious
                I_D_PT_spurious += cell_flux_PT_spurious; 
                I_D_DV_spurious += cell_flux_DV_spurious;
            }
                                            
            // Total 
            if (split == 1 || divTauOff == 1) { // tau_x term isn't accurated when splitting is on because cells are too stretched near the wall
                I_D_PT_tot += cell_flux_PT;
                I_D_DV_tot += cell_flux_DV;
            } else {
                I_D_PT_tot += cell_flux_PT;
                I_D_DV_tot += cell_flux_DV;
                I_D_DV_tot += (divTau_x * cellCenteredVariables[cell][CELL_AREA]) / V_inf; // with tau_x term correction
            }

            // Saving normalized (w.r.t. q_inf * c) flux value for each cell
            // Wave
            cellCenteredVariables[cell][CELL_FLUX_DRAG_PT_WV] = V_inf * cell_flux_PT_shock / (q_inf * chord);
            cellCenteredVariables[cell][CELL_FLUX_DRAG_DV_WV] = V_inf * cell_flux_DV_shock / (q_inf * chord);
            // Viscous
            cellCenteredVariables[cell][CELL_FLUX_DRAG_PT_VS] = V_inf * cell_flux_PT_viscous / (q_inf * chord);
            cellCenteredVariables[cell][CELL_FLUX_DRAG_DV_VS] = V_inf * cell_flux_DV_viscous / (q_inf * chord);
            // Spurious
            cellCenteredVariables[cell][CELL_FLUX_DRAG_PT_SP] = V_inf * cell_flux_PT_spurious / (q_inf * chord);
            cellCenteredVariables[cell][CELL_FLUX_DRAG_DV_SP] = V_inf * cell_flux_DV_spurious / (q_inf * chord);
            // Total value
            cellCenteredVariables[cell][CELL_FLUX_DRAG_PT] = V_inf * cell_flux_PT / (q_inf * chord);
            cellCenteredVariables[cell][CELL_FLUX_DRAG_DV] = V_inf * cell_flux_DV / (q_inf * chord);
        }

        // Compute irreversible drag 
        // Wave
        D_PT_wave = V_inf * I_D_PT_wave;
        D_DV_wave = V_inf * I_D_DV_wave;
        // Viscous
        D_PT_viscous = V_inf * I_D_PT_viscous;
        D_DV_viscous = V_inf * I_D_DV_viscous;
        // Spurious
        D_PT_spurious = V_inf * I_D_PT_spurious;
        D_DV_spurious = V_inf * I_D_DV_spurious;
        // Total
        D_PT_tot = V_inf * I_D_PT_tot;
        if (split == 1 || divTauOff == 1) { // tau_x term isn't accurated when splitting is on because cells are too stretched near the wall
            D_DV_tot = V_inf * I_D_DV_tot;
        } else {
            D_DV_tot = V_inf * I_D_DV_tot + Fx_f;
        }

        // Drag coefficients
        // Wave
        double Cd_wave_PT = D_PT_wave / (q_inf * chord);
        double Cd_wave_DV = D_DV_wave / (q_inf * chord);
        // Viscous
        double Cd_viscous_PT = D_PT_viscous / (q_inf * chord);
        double Cd_viscous_DV = D_DV_viscous / (q_inf * chord);
        // Spurious
        double Cd_spurious_PT = D_PT_spurious / (q_inf * chord);
        double Cd_spurious_DV = D_DV_spurious / (q_inf * chord);
        // Total
        double Cd_PT = D_PT_tot / (q_inf * chord);
        double Cd_DV = D_DV_tot / (q_inf * chord);   
        
        /* -------------------- LAMB VECTOR METHOD (THD-BASED APPROACH): F_ff = F_l + F_mrho + F_s + F_mu -------------------- */ 

        // Compute F_l components 
        for (int64_t cell = 0; cell < numCells; cell++) {

            double dF_l_x = 0.0; 
            double dF_l_y = 0.0;
            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];
            double dFl_std_x = 0.0;
            double dFl_std_y = 0.0;
    
            if (centroid_x > x_min && centroid_x < Xmax && centroid_y > y_min && centroid_y < y_max) {
                double cellArea = cellCenteredVariables[cell][CELL_AREA];
                if (vrt_std_on == 1) { // Standard Lamb Vector method (l = omega x V) (Wu et al.)
                    // Compute flow vorticity
                    double omega = cellCenteredVariables[cell][DV_DX] - cellCenteredVariables[cell][DU_DY]; 
                    // Calculate cell contribution
                    dFl_std_x = - omega * ccVariables[cell][V_CC] * ccVariables[cell][RHO_CC];
                    dFl_std_y =   omega * ccVariables[cell][U_CC] * ccVariables[cell][RHO_CC];
                    // Calculate volume integral
                    Fl_std_x += - dFl_std_x * cellArea;
                    Fl_std_y += - dFl_std_y * cellArea;
                }    
                // Unified approach (Minervino) --> F_l = volume integral of {rho*(T*gradS - gradH) + div(tau_v)}
                double rho_cc = ccVariables[cell][RHO_CC];
                double T_cc = ccVariables[cell][T_CC];
                double dSdx = cellCenteredVariables[cell][DS_DX];
                double dSdy = cellCenteredVariables[cell][DS_DY];
                double dHdx = cellCenteredVariables[cell][DH_DX];
                double dHdy = cellCenteredVariables[cell][DH_DY];
                double divTau_x = ccVariables[cell][DTAU11DX_FF] + ccVariables[cell][DTAU12DY_FF];
                double divTau_y = ccVariables[cell][DTAU12DX_FF] + ccVariables[cell][DTAU22DY_FF];
                
                // Calculate cell contribution
                dF_l_x = rho_cc * (T_cc * dSdx - dHdx) + divTau_x;
                dF_l_y = rho_cc * (T_cc * dSdy - dHdy) + divTau_y; 
                // Calculate volume integral
                F_l_x += - dF_l_x * cellArea;
                F_l_y += - dF_l_y * cellArea;
            }
        }

        // Compute F_mrho components: F_mrho = - volume integral of {r x [grad_rho x grad_K]}
        for (int cell = 0; cell < numCells; cell++) {

            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];

            if (centroid_x > x_min && centroid_x < Xmax && centroid_y > y_min && centroid_y < y_max) {

                // r vector components
                double r_x = ccVariables[cell][R_XC];
                double r_y = ccVariables[cell][R_YC];
                // Gradients from precomputed values
                double drho_dx = cellCenteredVariables[cell][DRHO_DX];
                double drho_dy = cellCenteredVariables[cell][DRHO_DY];
                double dK_dx = cellCenteredVariables[cell][DK_DX];
                double dK_dy = cellCenteredVariables[cell][DK_DY];
                double cellArea = cellCenteredVariables[cell][CELL_AREA];
                // Cross product A = grad_rho x grad_K
                double A = drho_dx * dK_dy - drho_dy * dK_dx;
                // Double cross product B = r x A (integrand)
                double B_x = A * r_y;
                double B_y = - A * r_x;
                // Volume integral contribution
                F_mrho_x += - B_x * cellArea;
                F_mrho_y += - B_y * cellArea;
            }
        }

        // Compute F_s components: 
        for (int64_t cell = 0; cell < numCells; cell++) {

            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];

            for (int localFace = 0; localFace < 4; localFace++) {
                // Face nodes index
                int64_t node1 = connectivity[cell][localFace];
                int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                if (node1 == node2) {
                    continue;
                }

                // Face centered interpolation of cell centered variables
                // Face centroid
                double x_f = 0.5 * (nodeVariables[node1 - 1][0] + nodeVariables[node2 - 1][0]);
                double y_f = 0.5 * (nodeVariables[node1 - 1][1] + nodeVariables[node2 - 1][1]);

                // Cell centered normals (no interpolation)
                double nx = normals[cell][localFace*2];
                double ny = normals[cell][localFace*2 + 1];

                // Looking for adjacent cell
                pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                int adjacentCell = -1; // none
                auto it = faceToCells.find(face); // reference to faceToCells
                if (it != faceToCells.end() && it->second.size() == 2) { // check if 'face' is shared exactly by two cells
                    adjacentCell = (it->second[0] == cell) ? it->second[1] : it->second[0]; // define adjacent cell (different from current cell)
                }

                double dsdx_f = 0.0, dsdy_f = 0.0, dHdx_f = 0.0, dHdy_f = 0.0, divTaux_f = 0.0, divTauy_f = 0.0; 
                double rholx_f = 0.0, rholy_f = 0.0, rholx_n = 0.0, rholx_p = 0.0, rholy_n = 0.0, rholy_p = 0.0;

                if(adjacentCell != -1) { // if adjacent cell exists
                    double dsdx_n = cellCenteredVariables[cell][DS_DX];
                    double dsdx_p = cellCenteredVariables[adjacentCell][DS_DX];
                    double dsdy_n = cellCenteredVariables[cell][DS_DY];
                    double dsdy_p = cellCenteredVariables[adjacentCell][DS_DY];
                    double dHdx_n = cellCenteredVariables[cell][DH_DX];
                    double dHdx_p = cellCenteredVariables[adjacentCell][DH_DX];
                    double dHdy_n = cellCenteredVariables[cell][DH_DY];
                    double dHdy_p = cellCenteredVariables[adjacentCell][DH_DY];
                    double divTaux_n = ccVariables[cell][DTAU11DX_FF] + ccVariables[cell][DTAU12DY_FF];
                    double divTaux_p = ccVariables[adjacentCell][DTAU11DX_FF] + ccVariables[adjacentCell][DTAU12DY_FF];
                    double divTauy_n = ccVariables[cell][DTAU12DX_FF] + ccVariables[cell][DTAU22DY_FF];
                    double divTauy_p = ccVariables[adjacentCell][DTAU12DX_FF] + ccVariables[adjacentCell][DTAU22DY_FF];
                    if (vrt_std_on == 1) {
                        rholx_n = - (cellCenteredVariables[cell][DV_DX] - cellCenteredVariables[cell][DU_DY]) * ccVariables[cell][V_CC] * ccVariables[cell][RHO_CC];
                        rholx_p = - (cellCenteredVariables[adjacentCell][DV_DX] - cellCenteredVariables[adjacentCell][DU_DY])  * ccVariables[adjacentCell][V_CC] * ccVariables[adjacentCell][RHO_CC];
                        rholy_n =   (cellCenteredVariables[cell][DV_DX] - cellCenteredVariables[cell][DU_DY])  * ccVariables[cell][U_CC] * ccVariables[cell][RHO_CC];
                        rholy_p =   (cellCenteredVariables[adjacentCell][DV_DX] - cellCenteredVariables[adjacentCell][DU_DY])  * ccVariables[adjacentCell][U_CC] * ccVariables[adjacentCell][RHO_CC];
                    }
                    
                    // Compute a 
                    double a = 0.0;
                    double dist_n = sqrt(pow(x_f - cellCenteredVariables[cell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[cell][CENTROID_Y], 2));
                    double dist_p = sqrt(pow(x_f - cellCenteredVariables[adjacentCell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[adjacentCell][CENTROID_Y], 2));
                    a = dist_n / (dist_p + dist_n);
                
                    dsdx_f = a * dsdx_p + (1.0 - a) * dsdx_n;
                    dsdy_f = a * dsdy_p + (1.0 - a) * dsdy_n;
                    dHdx_f = a * dHdx_p + (1.0 - a) * dHdx_n;
                    dHdy_f = a * dHdy_p + (1.0 - a) * dHdy_n;
                    divTaux_f = a * divTaux_p + (1.0 - a) * divTaux_n;
                    divTauy_f = a * divTauy_p + (1.0 - a) * divTauy_n;
                    if (vrt_std_on == 1) {
                        rholx_f = a * rholx_p + (1.0 - a) * rholx_n;
                        rholy_f = a * rholy_p + (1.0 - a) * rholy_n;
                    }
                    
                } else {
                    dsdx_f = cellCenteredVariables[cell][DS_DX];
                    dsdy_f = cellCenteredVariables[cell][DS_DY];
                    dHdx_f = cellCenteredVariables[cell][DH_DX];
                    dHdy_f = cellCenteredVariables[cell][DH_DY];
                    divTaux_f = ccVariables[cell][DTAU11DX_FF] + ccVariables[cell][DTAU12DY_FF];
                    divTauy_f = ccVariables[cell][DTAU12DX_FF] + ccVariables[cell][DTAU22DY_FF];
                    if (vrt_std_on == 1) {
                        rholx_f = - (cellCenteredVariables[cell][DV_DX] - cellCenteredVariables[cell][DU_DY]) * cellCenteredVariables[cell][momentumY_id];
                        rholy_f = (cellCenteredVariables[cell][DV_DX] - cellCenteredVariables[cell][DU_DY]) * cellCenteredVariables[cell][momentumX_id];
                    }
                }

                // Cell centered unitary normals
                double dS_f = sqrt(nx*nx + ny*ny);
                double nx_hat = nx / dS_f;
                double ny_hat = ny / dS_f;

                // Check if nodes are on the external boundaries of the control volume 
                if (nodeVariables[node1 - 1][numVars + 2] == 1 && nodeVariables[node2 - 1][numVars + 2] == 1) {
                    
                    // Compute face centered variables from nodal variables
                    // Face centered r vector components
                    double r_x1  = nodeVariables[node1 - 1][numVars + 4];
                    double r_y1  = nodeVariables[node1 - 1][numVars + 5];
                    double r_x2  = nodeVariables[node2 - 1][numVars + 4];
                    double r_y2  = nodeVariables[node2 - 1][numVars + 5];
                    double r_xf = 0.5 * (r_x1 + r_x2);
                    double r_yf = 0.5 * (r_y1 + r_y2);

                    // Face centered rho and T 
                    double rho1 = nodeVariables[node1 - 1][rho_id - 1];
                    double rho2 = nodeVariables[node2 - 1][rho_id - 1];
                    double rho_face = 0.5 * (rho1 + rho2);
                    double T1 = nodeVariables[node1 - 1][T_id - 1];
                    double T2 = nodeVariables[node2 - 1][T_id - 1];
                    double T_face = 0.5 * (T1 + T2);

                    if (centroid_x > x_min && centroid_x < Xmax && centroid_y > y_min && centroid_y < y_max) { // in this way we consider just the normals coming out of the external boundary
                        if (vrt_std_on == 1) { // Wu et al. formulation --> l = w x V
                            // Compute A = n x (rho*l)
                            double A = rholy_f * nx_hat - rholx_f * ny_hat;

                            // Compute r x A
                            double dFs_std_x = A * r_yf;
                            double dFs_std_y = - A * r_xf;

                            // Sum on faces and multiply for dS
                            Fs_std_x += - dFs_std_x * dS_f;
                            Fs_std_y += - dFs_std_y * dS_f;
                        }   
                        
                        // Minervino formulation --> F_s = - surface integral of {r x {n x [rho*(T*grad_S - grad_H) + div(tau_v)]}}dS
                        // Compute A = rho(T * gradS - gradH) + div(tau_v)
                        double A_x = rho_face * (T_face * dsdx_f - dHdx_f) + divTaux_f;
                        double A_y = rho_face * (T_face * dsdy_f - dHdy_f) + divTauy_f;

                        // Compute cross product B = n_hat x A
                        double B = nx_hat * A_y - ny_hat * A_x;
                        
                        // Compute cross product r x B
                        double dF_s_x = B * r_yf;
                        double dF_s_y = - B * r_xf;

                        // Sum on faces and multiply for dS
                        F_s_x += - dF_s_x * dS_f;
                        F_s_y += - dF_s_y * dS_f;
                    }
                }
            }
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Compute F_s as volume integral in order to obtain drag breakdown
        
        // 1) Compute gradients of rho(r dot l) and rho*r*l tensor using Green-Gauss method
        for (int64_t cell = 0; cell < numCells; cell++) {

            double grad_rdotl_std[2] = {0.0, 0.0};
            double grad_rdotl_thd[2] = {0.0, 0.0};
            double div_r_rhol_std[2] = {0.0, 0.0};
            double div_r_rhol_thd[2] = {0.0, 0.0};
            double cellArea = cellCenteredVariables[cell][CELL_AREA];
            
            for (int localFace = 0; localFace < 4; localFace++) {
                int64_t node1 = connectivity[cell][localFace];
                int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                if (node1 == node2) continue;

                double n_x = normals[cell][localFace * 2];
                double n_y = normals[cell][localFace * 2 + 1];

                // Face centroid
                double x_f = 0.5 * (nodeVariables[node1 - 1][0] + nodeVariables[node2 - 1][0]);
                double y_f = 0.5 * (nodeVariables[node1 - 1][1] + nodeVariables[node2 - 1][1]);

                // Looking for adjacent cell
                pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                int adjacentCell = -1; // none
                auto it = faceToCells.find(face); // reference to faceToCells
                if (it != faceToCells.end() && it->second.size() == 2) { // check if 'face' is shared exactly by two cells
                    adjacentCell = (it->second[0] == cell) ? it->second[1] : it->second[0]; // define adjacent cell (different from current cell)
                }

                // Face centered interpolation of cell centered variables
                double r_dot_l_std_f = 0.0, r_dot_l_thd_f = 0.0;
                double r_rhol_std_f[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
                double r_rhol_thd_f[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
                if (adjacentCell != -1) { // if adjacent cell exists
                    double r_dot_l_std_n = ccVariables[cell][R_DOT_L_std];
                    double r_dot_l_std_p = ccVariables[adjacentCell][R_DOT_L_std];
                    double r_dot_l_thd_n = ccVariables[cell][R_DOT_L_thd];
                    double r_dot_l_thd_p = ccVariables[adjacentCell][R_DOT_L_thd];
                    double r_rhol_std_n[2][2] = {
                        {ccVariables[cell][R_RHOL_std_XX], ccVariables[cell][R_RHOL_std_XY]},
                        {ccVariables[cell][R_RHOL_std_YX], ccVariables[cell][R_RHOL_std_YY]}
                    };
                    double r_rhol_std_p[2][2] = {
                        {ccVariables[adjacentCell][R_RHOL_std_XX], ccVariables[adjacentCell][R_RHOL_std_XY]},
                        {ccVariables[adjacentCell][R_RHOL_std_YX], ccVariables[adjacentCell][R_RHOL_std_YY]}
                    };
                    double r_rhol_thd_n[2][2] = {
                        {ccVariables[cell][R_RHOL_thd_XX], ccVariables[cell][R_RHOL_thd_XY]},
                        {ccVariables[cell][R_RHOL_thd_YX], ccVariables[cell][R_RHOL_thd_YY]}
                    };
                    double r_rhol_thd_p[2][2] = {
                        {ccVariables[adjacentCell][R_RHOL_thd_XX], ccVariables[adjacentCell][R_RHOL_thd_XY]},
                        {ccVariables[adjacentCell][R_RHOL_thd_YX], ccVariables[adjacentCell][R_RHOL_thd_YY]}
                    };

                    // Compute a
                    double a = 0.0;
                    double dist_n = sqrt(pow(x_f - cellCenteredVariables[cell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[cell][CENTROID_Y], 2));
                    double dist_p = sqrt(pow(x_f - cellCenteredVariables[adjacentCell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[adjacentCell][CENTROID_Y], 2));
                    a = dist_n / (dist_p + dist_n);

                    // Interpolate face centered variables
                    r_dot_l_std_f = a * r_dot_l_std_p + (1.0 - a) * r_dot_l_std_n;
                    r_dot_l_thd_f = a * r_dot_l_thd_p + (1.0 - a) * r_dot_l_thd_n;
                    r_rhol_std_f[0][0] = a * r_rhol_std_p[0][0] + (1.0 - a) * r_rhol_std_n[0][0];
                    r_rhol_std_f[0][1] = a * r_rhol_std_p[0][1] + (1.0 - a) * r_rhol_std_n[0][1];
                    r_rhol_std_f[1][0] = a * r_rhol_std_p[1][0] + (1.0 - a) * r_rhol_std_n[1][0];
                    r_rhol_std_f[1][1] = a * r_rhol_std_p[1][1] + (1.0 - a) * r_rhol_std_n[1][1];
                    r_rhol_thd_f[0][0] = a * r_rhol_thd_p[0][0] + (1.0 - a) * r_rhol_thd_n[0][0];
                    r_rhol_thd_f[0][1] = a * r_rhol_thd_p[0][1] + (1.0 - a) * r_rhol_thd_n[0][1];
                    r_rhol_thd_f[1][0] = a * r_rhol_thd_p[1][0] + (1.0 - a) * r_rhol_thd_n[1][0];
                    r_rhol_thd_f[1][1] = a * r_rhol_thd_p[1][1] + (1.0 - a) * r_rhol_thd_n[1][1];
                } else { // Do nothing as the lamb vector is 0 at the body wall!!
                    /* r_dot_l_std_f = ccVariables[cell][R_DOT_L_std];
                    r_dot_l_thd_f = ccVariables[cell][R_DOT_L_thd];
                    r_rhol_std_f[0][0] = ccVariables[cell][R_RHOL_std_XX];
                    r_rhol_std_f[0][1] = ccVariables[cell][R_RHOL_std_XY];
                    r_rhol_std_f[1][0] = ccVariables[cell][R_RHOL_std_YX];
                    r_rhol_std_f[1][1] = ccVariables[cell][R_RHOL_std_YY];
                    r_rhol_thd_f[0][0] = ccVariables[cell][R_RHOL_thd_XX];
                    r_rhol_thd_f[0][1] = ccVariables[cell][R_RHOL_thd_XY];
                    r_rhol_thd_f[1][0] = ccVariables[cell][R_RHOL_thd_YX];
                    r_rhol_thd_f[1][1] = ccVariables[cell][R_RHOL_thd_YY]; */
                }
                // Gradient components
                grad_rdotl_std[0] += r_dot_l_std_f * n_x;
                grad_rdotl_std[1] += r_dot_l_std_f * n_y;
                grad_rdotl_thd[0] += r_dot_l_thd_f * n_x;
                grad_rdotl_thd[1] += r_dot_l_thd_f * n_y;
                div_r_rhol_std[0] += (r_rhol_std_f[0][0] * n_x + r_rhol_std_f[1][0] * n_y);
                div_r_rhol_std[1] += (r_rhol_std_f[0][1] * n_x + r_rhol_std_f[1][1] * n_y);
                div_r_rhol_thd[0] += (r_rhol_thd_f[0][0] * n_x + r_rhol_thd_f[1][0] * n_y);
                div_r_rhol_thd[1] += (r_rhol_thd_f[0][1] * n_x + r_rhol_thd_f[1][1] * n_y);
            }

            // Normalize w.r.t cell area
            grad_rdotl_std[0] /= cellArea;
            grad_rdotl_std[1] /= cellArea;
            grad_rdotl_thd[0] /= cellArea;
            grad_rdotl_thd[1] /= cellArea;
            div_r_rhol_std[0] /= cellArea;
            div_r_rhol_std[1] /= cellArea;
            div_r_rhol_thd[0] /= cellArea;
            div_r_rhol_thd[1] /= cellArea;
            // Store computed values in auxiliary structures
            ccVariables[cell][GRAD_RDOTL_std_X] = grad_rdotl_std[0];
            ccVariables[cell][GRAD_RDOTL_std_Y] = grad_rdotl_std[1];
            ccVariables[cell][GRAD_RDOTL_thd_X] = grad_rdotl_thd[0];
            ccVariables[cell][GRAD_RDOTL_thd_Y] = grad_rdotl_thd[1];
            ccVariables[cell][DIV_R_RHOL_std_X] = div_r_rhol_std[0];
            ccVariables[cell][DIV_R_RHOL_std_Y] = div_r_rhol_std[1];
            ccVariables[cell][DIV_R_RHOL_thd_X] = div_r_rhol_thd[0];
            ccVariables[cell][DIV_R_RHOL_thd_Y] = div_r_rhol_thd[1];
        }

        // 2) Compute F_s as volume integral of (- grad_rho(r dot l)) + div(rho*r*l) + breakdown of F_s contribution to parasite drag 
        double F_s_x_std = 0.0, F_s_y_std = 0.0, F_s_x_thd = 0.0, F_s_y_thd = 0.0; // volume integral
        double F_s_x_std_visc = 0.0, F_s_x_std_wave = 0.0, F_s_x_std_spu = 0.0;
        double F_s_x_thd_visc = 0.0, F_s_x_thd_wave = 0.0, F_s_x_thd_spu = 0.0;

        for(int64_t cell = 0; cell < numCells; cell++) {
            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];            
            double cellArea = cellCenteredVariables[cell][CELL_AREA];

            if (centroid_x > x_min && centroid_x < Xmax && centroid_y > y_min && centroid_y < y_max) {

                double grad_rdotl_std_x = ccVariables[cell][GRAD_RDOTL_std_X];
                double grad_rdotl_thd_x = ccVariables[cell][GRAD_RDOTL_thd_X];
                double div_r_rhol_std_x = ccVariables[cell][DIV_R_RHOL_std_X];
                double div_r_rhol_thd_x = ccVariables[cell][DIV_R_RHOL_thd_X]; 

                double term_std = (- grad_rdotl_std_x + div_r_rhol_std_x);
                double term_thd = (- grad_rdotl_thd_x + div_r_rhol_thd_x);

                // Calculate viscous, wave and spurious F_s contributions
                if(cellCenteredVariables[cell][CELL_TYPE] == 1.0){ // Compute wave contribution
                    F_s_x_std_wave += term_std * cellArea;
                    F_s_x_thd_wave += term_thd * cellArea;
                } else if(cellCenteredVariables[cell][CELL_TYPE] == 2.0){ // Compute viscous contribution
                    F_s_x_std_visc += term_std * cellArea;
                    F_s_x_thd_visc += term_thd * cellArea;
                } else if(cellCenteredVariables[cell][CELL_TYPE] == 0.0){ // Compute spurious contribution
                    F_s_x_std_spu += term_std * cellArea;
                    F_s_x_thd_spu += term_thd * cellArea;
                }
                // Total contribution
                F_s_x_std += term_std * cellArea;
                F_s_x_thd += term_thd * cellArea;
                // Normalize and save volumetric cell contribution to display in Tecplot parasite drag contributions, i.e F_s 
                if(cellCenteredVariables[cell][CELL_TYPE] == 1.0){ // Wave
                    cellCenteredVariables[cell][FS_STD_WV] = term_std * cellArea / (q_inf * chord);
                    cellCenteredVariables[cell][FS_THD_WV] = term_thd * cellArea / (q_inf * chord);
                } else if(cellCenteredVariables[cell][CELL_TYPE] == 2.0){ // Viscous
                    cellCenteredVariables[cell][FS_STD_VS] = term_std * cellArea / (q_inf * chord);
                    cellCenteredVariables[cell][FS_THD_VS] = term_thd * cellArea / (q_inf * chord);
                } else if(cellCenteredVariables[cell][CELL_TYPE] == 0.0){ // Spurious
                    cellCenteredVariables[cell][FS_STD_SP] = term_std * cellArea / (q_inf * chord);
                    cellCenteredVariables[cell][FS_THD_SP] = term_thd * cellArea / (q_inf * chord);
                }
                // Total value
                cellCenteredVariables[cell][FS_STD] = term_std * cellArea / (q_inf * chord);
                cellCenteredVariables[cell][FS_THD] = term_thd * cellArea/ (q_inf * chord);
            }
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        // F_mud THD components: F_mu = surface integral of {r x [n x div(tau)]dS} + surface integral of {[tau_v dot n]dS} 
        double dFmu1_x = 0.0, dFmu2_x = 0.0, dFmu1_y = 0.0, dFmu2_y = 0.0;
        for (int64_t cell = 0; cell < numCells; cell++) { 

            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];

            for (int localFace = 0; localFace < 4; localFace++) {
                // Face nodes index
                int64_t node1 = connectivity[cell][localFace];
                int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                if (node1 == node2) {
                    continue;
                }

                // Cell centered unitary normals
                double nx = normals[cell][localFace * 2];
                double ny = normals[cell][localFace * 2 + 1];

                // Looking for adjacent cell
                pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                int adjacentCell = -1; // none
                auto it = faceToCells.find(face); // reference to faceToCells
                if (it != faceToCells.end() && it->second.size() == 2) { // check if 'face' is shared exactly by two cells
                    adjacentCell = (it->second[0] == cell) ? it->second[1] : it->second[0]; // define adjacent cell (different from current cell)
                }

                // Face centered interpolation of cell variables
                double tau11_f = 0.0, tau12_f = 0.0, tau22_f = 0.0, divTaux_f = 0.0, divTauy_f = 0.0;

                if(adjacentCell != -1) { // if adjacent cell exists
                    double tau11_n = cellCenteredVariables[cell][TAU11_FF];
                    double tau11_p = cellCenteredVariables[adjacentCell][TAU11_FF];
                    double tau12_n = cellCenteredVariables[cell][TAU12_FF];
                    double tau12_p = cellCenteredVariables[adjacentCell][TAU12_FF];
                    double tau22_n = cellCenteredVariables[cell][TAU22_FF];
                    double tau22_p = cellCenteredVariables[adjacentCell][TAU22_FF];
                    double divTaux_n = ccVariables[cell][DTAU11DX_FF] + ccVariables[cell][DTAU12DY_FF];
                    double divTaux_p = ccVariables[adjacentCell][DTAU11DX_FF] + ccVariables[adjacentCell][DTAU12DY_FF];
                    double divTauy_n = ccVariables[cell][DTAU12DX_FF] + ccVariables[cell][DTAU22DY_FF];
                    double divTauy_p = ccVariables[adjacentCell][DTAU12DX_FF] + ccVariables[adjacentCell][DTAU22DY_FF];

                    // Face centroid
                    double x_f = 0.5 * (nodeVariables[node1 - 1][0] + nodeVariables[node2 - 1][0]);
                    double y_f = 0.5 * (nodeVariables[node1 - 1][1] + nodeVariables[node2 - 1][1]);

                    // Compute a 
                    double a = 0.0;
                    double dist_n = sqrt(pow(x_f - cellCenteredVariables[cell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[cell][CENTROID_Y], 2));
                    double dist_p = sqrt(pow(x_f - cellCenteredVariables[adjacentCell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[adjacentCell][CENTROID_Y], 2));
                    a = dist_n / (dist_p + dist_n);
                
                    tau11_f = a * tau11_p + (1.0 - a) * tau11_n;
                    tau12_f = a * tau12_p + (1.0 - a) * tau12_n;
                    tau22_f = a * tau22_p + (1.0 - a) * tau22_n;
                    divTaux_f = a * divTaux_p + (1.0 - a) * divTaux_n;
                    divTauy_f = a * divTauy_p + (1.0 - a) * divTauy_n;
                } else {
                    tau11_f= cellCenteredVariables[cell][TAU11_FF];
                    tau12_f = cellCenteredVariables[cell][TAU12_FF];
                    tau22_f = cellCenteredVariables[cell][TAU22_FF];
                    divTaux_f = ccVariables[cell][DTAU11DX_FF] + ccVariables[cell][DTAU12DY_FF];
                    divTauy_f = ccVariables[cell][DTAU12DX_FF] + ccVariables[cell][DTAU22DY_FF];
                }

                // Cell centered unitary normals
                double dS_f = sqrt(nx*nx + ny*ny);
                double nx_hat = nx / dS_f;
                double ny_hat = ny / dS_f; 

                // Check if nodes are on the external boundary of the control volume
                if (nodeVariables[node1 - 1][numVars + 2] == 1 && nodeVariables[node2 - 1][numVars + 2] == 1) {

                    // Face centered r vector components
                    double r_x1  = nodeVariables[node1 - 1][numVars + 4];
                    double r_y1  = nodeVariables[node1 - 1][numVars + 5];
                    double r_x2  = nodeVariables[node2 - 1][numVars + 4];
                    double r_y2  = nodeVariables[node2 - 1][numVars + 5];
                    double r_xf = 0.5 * (r_x1 + r_x2);
                    double r_yf = 0.5 * (r_y1 + r_y2);

                    if (centroid_x > x_min && centroid_x < Xmax && centroid_y > y_min && centroid_y < y_max) {

                        // Cross-product [n_hat x (div(tau))]
                        double A = nx_hat * divTauy_f - ny_hat * divTaux_f;

                        // Compute first integral term (r x A)
                        dFmu1_x = A * r_yf;
                        dFmu1_y = - A * r_xf;

                        // Compute second integral term (tau_v*n)                                 
                        dFmu2_x = tau11_f * nx_hat + tau12_f * ny_hat;
                        dFmu2_y = tau12_f * nx_hat + tau22_f * ny_hat;

                        // Sum on faces and multiply for dS
                        F_mud_x_thd += (dFmu1_x + dFmu2_x) * dS_f;
                        F_mud_y_thd += (dFmu1_y + dFmu2_y) * dS_f;
                    }
                }
            }   
        } 

        // F_mu STD (with tau = tau_v + tau_R) components: F_mu = surface integral of {r x [n x div(tau)]dS} + surface integral of {[tau_v dot n]dS} 
        double dFmu1_x_std = 0.0, dFmu2_x_std = 0.0, dFmu1_y_std = 0.0, dFmu2_y_std = 0.0;
        for (int64_t cell = 0; cell < numCells; cell++) { 

            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];

            for (int localFace = 0; localFace < 4; localFace++) {
                // Face nodes index
                int64_t node1 = connectivity[cell][localFace];
                int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                if (node1 == node2) {
                    continue;
                }

                // Cell centered unitary normals
                double nx = normals[cell][localFace * 2];
                double ny = normals[cell][localFace * 2 + 1];

                // Looking for adjacent cell
                pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                int adjacentCell = -1; // none
                auto it = faceToCells.find(face); // reference to faceToCells
                if (it != faceToCells.end() && it->second.size() == 2) { // check if 'face' is shared exactly by two cells
                    adjacentCell = (it->second[0] == cell) ? it->second[1] : it->second[0]; // define adjacent cell (different from current cell)
                }

                // Face centered interpolation of cell variables
                double tau11_f = 0.0, tau12_f = 0.0, tau22_f = 0.0, divTaux_f = 0.0, divTauy_f = 0.0;

                if(adjacentCell != -1) { // if adjacent cell exists
                    double tau11_n = cellCenteredVariables[cell][TAU11_NF];
                    double tau11_p = cellCenteredVariables[adjacentCell][TAU11_NF];
                    double tau12_n = cellCenteredVariables[cell][TAU12_NF];
                    double tau12_p = cellCenteredVariables[adjacentCell][TAU12_NF];
                    double tau22_n = cellCenteredVariables[cell][TAU22_NF];
                    double tau22_p = cellCenteredVariables[adjacentCell][TAU22_NF];
                    double divTaux_n = ccVariables[cell][DTAU11DX_NF] + ccVariables[cell][DTAU12DY_NF];
                    double divTaux_p = ccVariables[adjacentCell][DTAU11DX_NF] + ccVariables[adjacentCell][DTAU12DY_NF];
                    double divTauy_n = ccVariables[cell][DTAU12DX_NF] + ccVariables[cell][DTAU22DY_NF];
                    double divTauy_p = ccVariables[adjacentCell][DTAU12DX_NF] + ccVariables[adjacentCell][DTAU22DY_NF];

                    // Face centroid
                    double x_f = 0.5 * (nodeVariables[node1 - 1][0] + nodeVariables[node2 - 1][0]);
                    double y_f = 0.5 * (nodeVariables[node1 - 1][1] + nodeVariables[node2 - 1][1]);

                    // Compute a 
                    double a = 0.0;
                    double dist_n = sqrt(pow(x_f - cellCenteredVariables[cell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[cell][CENTROID_Y], 2));
                    double dist_p = sqrt(pow(x_f - cellCenteredVariables[adjacentCell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[adjacentCell][CENTROID_Y], 2));
                    a = dist_n / (dist_p + dist_n);
                
                    tau11_f = a * tau11_p + (1.0 - a) * tau11_n;
                    tau12_f = a * tau12_p + (1.0 - a) * tau12_n;
                    tau22_f = a * tau22_p + (1.0 - a) * tau22_n;
                    divTaux_f = a * divTaux_p + (1.0 - a) * divTaux_n;
                    divTauy_f = a * divTauy_p + (1.0 - a) * divTauy_n;
                } else {
                    tau11_f= cellCenteredVariables[cell][TAU11_NF];
                    tau12_f = cellCenteredVariables[cell][TAU12_NF];
                    tau22_f = cellCenteredVariables[cell][TAU22_NF];
                    divTaux_f = ccVariables[cell][DTAU11DX_NF] + ccVariables[cell][DTAU12DY_NF];
                    divTauy_f = ccVariables[cell][DTAU12DX_NF] + ccVariables[cell][DTAU22DY_NF];
                }

                // Cell centered unitary normals
                double dS_f = sqrt(nx*nx + ny*ny);
                double nx_hat = nx / dS_f;
                double ny_hat = ny / dS_f; 

                // Check if nodes are on the external boundary of the control volume
                if (nodeVariables[node1 - 1][numVars + 2] == 1 && nodeVariables[node2 - 1][numVars + 2] == 1) {

                    // Face centered r vector components
                    double r_x1  = nodeVariables[node1 - 1][numVars + 4];
                    double r_y1  = nodeVariables[node1 - 1][numVars + 5];
                    double r_x2  = nodeVariables[node2 - 1][numVars + 4];
                    double r_y2  = nodeVariables[node2 - 1][numVars + 5];
                    double r_xf = 0.5 * (r_x1 + r_x2);
                    double r_yf = 0.5 * (r_y1 + r_y2);

                    if (centroid_x > x_min && centroid_x < Xmax && centroid_y > y_min && centroid_y < y_max) {

                        // Cross-product [n_hat x (div(tau))]
                        double A = nx_hat * divTauy_f - ny_hat * divTaux_f;

                        // Compute first integral term (r x A)
                        dFmu1_x_std = A * r_yf;
                        dFmu1_y_std = - A * r_xf;

                        // Compute second integral term (tau_v*n)                                 
                        dFmu2_x_std = tau11_f * nx_hat + tau12_f * ny_hat;
                        dFmu2_y_std = tau12_f * nx_hat + tau22_f * ny_hat;

                        // Sum on faces and multiply for dS
                        F_mu_x_std += (dFmu1_x_std + dFmu2_x_std) * dS_f;
                        F_mu_y_std += (dFmu1_y_std + dFmu2_y_std) * dS_f;
                    }
                }
            }   
        }

        /* // F_mut turbulent correction term: - volume integral of {r x [nabla x (rho*grad_Kt)]}
        for (int cell = 0; cell < numCells; cell++) {

            double rotRhoDKt[2] = {0.0, 0.0};
            double cellArea = cellCenteredVariables[cell][CELL_AREA];
            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];

            for (int localFace = 0; localFace < 4; localFace++) {
                int64_t node1 = connectivity[cell][localFace];
                int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                if (node1 == node2) continue; 

                // Normals
                double n_x = normals[cell][localFace * 2];
                double n_y = normals[cell][localFace * 2 + 1];

                // Face centroid
                double x_f = 0.5 * (nodeVariables[node1 - 1][0] + nodeVariables[node2 - 1][0]);
                double y_f = 0.5 * (nodeVariables[node1 - 1][1] + nodeVariables[node2 - 1][1]);

                // Looking for adjacent cell
                pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                int adjacentCell = -1; // none
                auto it = faceToCells.find(face); // reference to faceToCells
                if (it != faceToCells.end() && it->second.size() == 2) { // check if 'face' is shared exactly by two cells
                    adjacentCell = (it->second[0] == cell) ? it->second[1] : it->second[0]; // define adjacent cell (different from current cell)
                }

                // Face centered values of rho and dKtdx,dKtdy for current cell (n) and adjacent cell (p)
                double rho_face = 0.0, dKtdx_face = 0.0, dKtdy_face = 0.0, rho_dKtdx_face = 0.0, rho_dKtdy_face = 0.0;
                rho_face = 0.5 * (nodeVariables[node1 - 1][rho_id - 1] + nodeVariables[node2 - 1][rho_id - 1]) ;

                // Face centered interpolation
                if(adjacentCell != -1) { // if adjacent cell exists

                    double dKtdx_n = cellCenteredVariables[cell][DKT_DX];
                    double dKtdy_n = cellCenteredVariables[cell][DKT_DY];
                    double dKtdx_p = cellCenteredVariables[adjacentCell][DKT_DX];
                    double dKtdy_p = cellCenteredVariables[adjacentCell][DKT_DY];

                    // Compute a 
                    double a = 0.0;
                    double dist_n = sqrt(pow(x_f - cellCenteredVariables[cell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[cell][CENTROID_Y], 2));
                    double dist_p = sqrt(pow(x_f - cellCenteredVariables[adjacentCell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[adjacentCell][CENTROID_Y], 2));
                    a = dist_n / (dist_p + dist_n);
                    
                    dKtdx_face = a * dKtdx_p + (1.0 - a) * dKtdx_n;
                    dKtdy_face = a * dKtdy_p + (1.0 - a) * dKtdy_n;

                    rho_dKtdx_face = rho_face * dKtdx_face;
                    rho_dKtdy_face = rho_face * dKtdy_face;

                } else {

                    double dKtdx_n = cellCenteredVariables[cell][DKT_DX];
                    double dKtdy_n = cellCenteredVariables[cell][DKT_DY];
                    double rhodKtdx_n = rho_face * dKtdx_n;
                    double rhodKtdy_n = rho_face * dKtdy_n;

                    rho_dKtdx_face = rhodKtdx_n;
                    rho_dKtdy_face = rhodKtdy_n;
                }

                // Compute cell centered values of (rho * dKtdx) and (rho * dKtdy)
                rotRhoDKt[0] += rho_dKtdx_face * n_y; // ddy(rho * dKtdx)
                rotRhoDKt[1] += rho_dKtdy_face * n_x; // ddx(rho * dKtdy)
            }

            // Normalize w.r.t cell area
            rotRhoDKt[0] /= cellArea; // ddy(rho * dKtdx)
            rotRhoDKt[1] /= cellArea; // ddx(rho * dKtdy)

            cellCenteredVariables[cell][DDY_RHODKTDX] = rotRhoDKt[0];
            cellCenteredVariables[cell][DDX_RHODKTDY] = rotRhoDKt[1];
        }

        // Compute F_mut volume integral 
        for (int cell = 0; cell < numCells; cell++) {
            
            double B_x = 0.0;
            double B_y = 0.0;
            double A = 0.0;
            double centroid_x = cellCenteredVariables[cell][centroidX_id];
            double centroid_y = cellCenteredVariables[cell][centroidY_id];

            if (centroid_x > x_min && centroid_x < Xmax && centroid_y > y_min && centroid_y < y_max) {

                // r vector components
                double r_x = ccVariables[cell][R_XC];
                double r_y = ccVariables[cell][R_YC];
                // Cell area
                double cellArea = cellCenteredVariables[cell][CELL_AREA];

                // Cross product A = nabla x (rho * gradKt) = rot(rho*gradKt)
                A = cellCenteredVariables[cell][DDX_RHODKTDY] - cellCenteredVariables[cell][DDY_RHODKTDX];

                // Double cross product B = r x A (integrand)
                B_x = A * r_y;
                B_y = - A * r_x;

                // Volume integral contribution
                F_mut_x += - B_x * cellArea;
                F_mut_y += - B_y * cellArea;
            }
        } */

        // F_mut --> DMTs: F_mut = volume integral of (rho*gradKt) - surface integral of (r x n x (rho*gradKt))
        // 1) Volume integral
        for (int cell = 0; cell < numCells; cell++) {

            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];
            double rhoDkt_x = 0.0, rhoDkt_y = 0.0;

            if (centroid_x > x_min && centroid_x < Xmax && centroid_y > y_min && centroid_y < y_max) {

                // Cell centered values
                double rho_cc = ccVariables[cell][RHO_CC];
                double dkt_dx = cellCenteredVariables[cell][DKT_DX];
                double dkt_dy = cellCenteredVariables[cell][DKT_DY];
                double cellArea = cellCenteredVariables[cell][CELL_AREA];
                rhoDkt_x = rho_cc * dkt_dx;
                rhoDkt_y = rho_cc * dkt_dy;
                // Volume integral contribution
                F_mut1_x += - rhoDkt_x * cellArea;
                F_mut1_y += - rhoDkt_y * cellArea;
            }
        }

        // 2) Surface integral
        for (int64_t cell = 0; cell < numCells; cell++) { 

            double centroid_x = cellCenteredVariables[cell][CENTROID_X];
            double centroid_y = cellCenteredVariables[cell][CENTROID_Y];

            for (int localFace = 0; localFace < 4; localFace++) {
                int64_t node1 = connectivity[cell][localFace];
                int64_t node2 = connectivity[cell][(localFace + 1) % 4];

                if (node1 == node2) {
                    continue;
                }  

                // Normals
                double nx = normals[cell][localFace * 2];
                double ny = normals[cell][localFace * 2 + 1];

                // Face centroid
                double x_f = 0.5 * (nodeVariables[node1 - 1][0] + nodeVariables[node2 - 1][0]);
                double y_f = 0.5 * (nodeVariables[node1 - 1][1] + nodeVariables[node2 - 1][1]);

                // Looking for adjacent cell
                pair<int64_t, int64_t> face = {min(node1, node2), max(node1, node2)};
                int adjacentCell = -1; // none
                auto it = faceToCells.find(face); // reference to faceToCells
                if (it != faceToCells.end() && it->second.size() == 2) { // check if 'face' is shared exactly by two cells
                    adjacentCell = (it->second[0] == cell) ? it->second[1] : it->second[0]; // define adjacent cell (different from current cell)
                }

                // Face centered values of dKtdx,dKtdy for current cell (n) and adjacent cell (p)
                double dKtdx_face = 0.0, dKtdy_face = 0.0;
                if(adjacentCell != -1) { // if adjacent cell exists

                    double dKtdx_n = cellCenteredVariables[cell][DKT_DX];
                    double dKtdy_n = cellCenteredVariables[cell][DKT_DY];
                    double dKtdx_p = cellCenteredVariables[adjacentCell][DKT_DX];
                    double dKtdy_p = cellCenteredVariables[adjacentCell][DKT_DY];

                    // Compute a 
                    double a = 0.0;
                    double dist_n = sqrt(pow(x_f - cellCenteredVariables[cell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[cell][CENTROID_Y], 2));
                    double dist_p = sqrt(pow(x_f - cellCenteredVariables[adjacentCell][CENTROID_X], 2) + pow(y_f - cellCenteredVariables[adjacentCell][CENTROID_Y], 2));
                    a = dist_n / (dist_p + dist_n);
                    
                    dKtdx_face = a * dKtdx_p + (1.0 - a) * dKtdx_n;
                    dKtdy_face = a * dKtdy_p + (1.0 - a) * dKtdy_n;
                } else {
                    dKtdx_face = cellCenteredVariables[cell][DKT_DX];
                    dKtdy_face = cellCenteredVariables[cell][DKT_DY];
                }

                // Cell centered unitary normals
                double dS_f = sqrt(nx*nx + ny*ny);
                double nx_hat = nx / dS_f;
                double ny_hat = ny / dS_f; 

                // Check if nodes are on the external boundary of the control volume
                if (nodeVariables[node1 - 1][numVars + 2] == 1 && nodeVariables[node2 - 1][numVars + 2] == 1) {

                    // Face centered r vector components
                    double r_x1  = nodeVariables[node1 - 1][numVars + 4];
                    double r_y1  = nodeVariables[node1 - 1][numVars + 5];
                    double r_x2  = nodeVariables[node2 - 1][numVars + 4];
                    double r_y2  = nodeVariables[node2 - 1][numVars + 5];
                    double r_xf = 0.5 * (r_x1 + r_x2);
                    double r_yf = 0.5 * (r_y1 + r_y2);

                    // Face centered rho
                    double rho1 = nodeVariables[node1 - 1][rho_id - 1];
                    double rho2 = nodeVariables[node2 - 1][rho_id - 1];
                    double rho_face = 0.5 * (rho1 + rho2);

                    if (centroid_x > x_min && centroid_x < Xmax && centroid_y > y_min && centroid_y < y_max) {

                        // Compute A = rho * gradKt
                        double A_x = rho_face * dKtdx_face;
                        double A_y = rho_face * dKtdy_face;
                        // Compute cross product B = n_hat x A
                        double B = nx_hat * A_y - ny_hat * A_x;

                        // Compute cross product r x B
                        double dF_mut2_x = B * r_yf;
                        double dF_mut2_y = - B * r_xf;

                        // Sum on faces and multiply for dS
                        F_mut2_x += - dF_mut2_x * dS_f;
                        F_mut2_y += - dF_mut2_y * dS_f;
                    }
                }
            }   
        }

        // THD APPROACH - Total viscous contribution F_mu = F_mud + F_mut (laminar + turbulent correction)
        F_mu_x_thd = F_mud_x_thd + F_mut1_x + F_mut2_x;
        F_mu_y_thd = F_mud_y_thd + F_mut1_y + F_mut2_y;
        // STD APPROACH - Total viscous contribution F_mu (tau = tau_v + tau_R)
        // Already computed

        // Compute Cd anc Cl contributions
        double F_ff_y = 0.0, F_ff_x = 0.0, Cl_vf_std = 0.0, Cd_vf_std = 0.0, Cl_vf_thd = 0.0, Cd_vf_thd = 0.0;
        double Cl_compr = 0.0, Cd_compr = 0.0, Cl_fs_std = 0.0, Cd_fs_std = 0.0, Cl_fs_thd = 0.0, Cd_fs_thd = 0.0;
        double Cl_fmu_thd = 0.0, Cd_fmu_thd = 0.0, Cl_fmu_std = 0.0, Cd_fmu_std = 0.0, Cl_vrt_thd = 0.0, Cd_vrt_thd = 0.0, Cl_vrt_std = 0.0, Cd_vrt_std= 0.0;
        // Minervino unified thd approach
        Cl_vf_thd = F_l_y / (q_inf * chord);
        Cd_vf_thd = F_l_x / (q_inf * chord);
        Cl_compr = F_mrho_y / (q_inf * chord);
        Cd_compr = F_mrho_x / (q_inf * chord);
        Cl_fs_thd = F_s_y / (q_inf * chord);
        Cd_fs_thd = F_s_x / (q_inf * chord);
        Cl_fmu_thd = F_mu_y_thd / (q_inf * chord);
        Cd_fmu_thd = F_mu_x_thd / (q_inf * chord);
        F_ff_y = F_l_y + F_mrho_y + F_s_y + F_mu_y_thd;
        F_ff_x = F_l_x + F_mrho_x + F_s_x + F_mu_x_thd;
        Cl_vrt_thd = F_ff_y / (q_inf * chord);
        Cd_vrt_thd = F_ff_x / (q_inf * chord);
        // Wu et al. standard formulation
        if (vrt_std_on == 1) {
            Cl_vf_std = Fl_std_y / (q_inf * chord);
            Cd_vf_std = Fl_std_x / (q_inf * chord);
            Cl_fs_std = Fs_std_y / (q_inf * chord);
            Cd_fs_std = Fs_std_x / (q_inf * chord);
            Cl_fmu_std = F_mu_y_std / (q_inf * chord);
            Cd_fmu_std = F_mu_x_std / (q_inf * chord);
            Cl_vrt_std = (Fl_std_y + F_mrho_y + Fs_std_y + F_mu_y_std) / (q_inf * chord);
            Cd_vrt_std = (Fl_std_x + F_mrho_x + Fs_std_x + F_mu_x_std) / (q_inf * chord);
        }

        // Calculate induced drag 
        double Cd_i_std = 0.0, Cd_i_thd = 0.0;
        Cd_i_std = Cd_vf_std + Cd_compr;
        Cd_i_thd = Cd_vf_thd + Cd_compr;
        // Calculate parasite drag and perform drag breakdown
        double Cd_p_std = 0.0, Cd_p_thd = 0.0, Cd_wv_std = 0.0, Cd_wv_thd = 0.0, Cd_vs_std = 0.0, Cd_vs_thd = 0.0, Cd_sp_std = 0.0, Cd_sp_thd = 0.0;
        Cd_p_std = Cd_fmu_std + (F_s_x_std / (q_inf * chord)); 
        Cd_p_thd = Cd_fmu_thd + (F_s_x_thd / (q_inf * chord)); 
        Cd_wv_std = F_s_x_std_wave / (q_inf * chord);
        Cd_wv_thd = F_s_x_thd_wave / (q_inf * chord);
        Cd_vs_std = Cd_fmu_std + (F_s_x_std_visc / (q_inf * chord));
        Cd_vs_thd = Cd_fmu_thd + (F_s_x_thd_visc / (q_inf * chord));
        Cd_sp_std = F_s_x_std_spu / (q_inf * chord);
        Cd_sp_thd = F_s_x_thd_spu / (q_inf * chord);

        // Print drag breakdown
        std::cout << "\n" << endl;
        printf("%s\n", "|                                     FAR-FIELD AERODYNAMIC DRAG BREAKDOWN                                    |");
        printf("%s\n", "|-------------------------------------------------------------------------------------------------------------|");
        printf("%s%8.1f%s\n", "|                                          X_w/c: ", Xmax, "                                                    |");
        printf("%s\n", "|-------------------------------------------------------------------------------------------------------------|");
        
        // THD METHODS
        printf("%s\n", "|                                   THERMODYNAMIC METHODS                                                     |");
        printf("%s\n", "|-------------------------------------------------------------------------------------------------------------|");
        printf("%s%8.1f%s\n", "|  TOTAL DRAG COUNTS   (DESTARAC & VAN DER VOOREN)             : ", 10000 * Cd_DV, "                                     |");
        printf("%s%8.1f%s\n", "|                                ----------------------------- WAVE     : ", 10000 * Cd_wave_DV, "                            |");
        printf("%s%8.1f%s\n", "|                                -------------------------- VISCOUS     : ", 10000 * Cd_viscous_DV, "                            |");
        printf("%s%8.1f%s\n", "|                                ------------------------- SPURIOUS     : ", 10000 * Cd_spurious_DV, "                            |");
        printf("%s%8.1f%s\n", "|  TOTAL DRAG COUNTS   (PAPARONE & TOGNACCINI)                 : ", 10000 * Cd_PT, "                                     |");
        printf("%s%8.1f%s\n", "|                                ----------------------------- WAVE     : ", 10000 * Cd_wave_PT, "                            |");
        printf("%s%8.1f%s\n", "|                                -------------------------- VISCOUS     : ", 10000 * Cd_viscous_PT, "                            |");
        printf("%s%8.1f%s\n", "|                                ------------------------- SPURIOUS     : ", 10000 * Cd_spurious_PT, "                            |");
        printf("%s\n", "|-------------------------------------------------------------------------------------------------------------|");
        
        // VRT MINERVINO AND TOGNACCINI
        printf("%s\n", "|                         LAMB VECTOR METHOD - THD-BASED APPROACH (MINERVINO & TOGNACCCINI)                   |");
        printf("%s\n", "|-------------------------------------------------------------------------------------------------------------|");
        printf("%s%8.1f%s\n", "|  PARASITE DRAG                                               : ", 10000 * Cd_p_thd, "                                     |");
        printf("%s%8.1f%s\n", "|                                ----------------------------- WAVE     : ", 10000 * Cd_wv_thd, "                            |");
        printf("%s%8.1f%s\n", "|                                -------------------------- VISCOUS     : ", 10000 * Cd_vs_thd, "                            |");
        printf("%s%8.1f%s\n", "|                                ------------------------- SPURIOUS     : ", 10000 * Cd_sp_thd, "                            |");
        printf("%s%8.1f%s\n", "|                                -------------------------- INDUCED     : ", 10000 * Cd_i_thd, "                            |");
        printf("%s\n", "|-------------------------------------------------------------------------------------------------------------|");
        
        // VRT WU ET AL.
        if (vrt_std_on == 1) {
            printf("%s\n", "|                         LAMB VECTOR METHOD - STD APPROACH (WU et al.)                                       |");
            printf("%s\n", "|-------------------------------------------------------------------------------------------------------------|");
            printf("%s%8.1f%s\n", "|  PARASITE DRAG                                               : ", 10000 * Cd_p_std, "                                     |");
            printf("%s%8.1f%s\n", "|                                ----------------------------- WAVE     : ", 10000 * Cd_wv_std, "                            |");
            printf("%s%8.1f%s\n", "|                                -------------------------- VISCOUS     : ", 10000 * Cd_vs_std, "                            |");
            printf("%s%8.1f%s\n", "|                                ------------------------- SPURIOUS     : ", 10000 * Cd_sp_std, "                            |");
            printf("%s%8.1f%s\n", "|                                -------------------------- INDUCED     : ", 10000 * Cd_i_std, "                            |");
            printf("%s\n", "|-------------------------------------------------------------------------------------------------------------|");
        }

        // Print aerodynamic coefficients
        std::cout << "\n" << endl;
        printf("%s\n",					"|        FAR-FIELD AERODYNAMIC FORCE BREAK-DOWN  (LAMB VECTOR METHOD: MINERVINO & TOGNACCINI APPROACH)        |");
        printf("%s\n",					"|------------------------------------------------------------------------------------------------------------ |");
        printf("%s%8.1f%s\n",	        "|                                                 X_w/c: ", (Xmax), "\t\t\t\t\t"                          "      |");
        printf("%s\n",					"|-------------------------------------------------------------------------------------------------------------|");
        printf("%s\n",					"|                                                                                                             |");
        printf("%s%8.1f%s\n",			"|  TOTAL DRAG COUNTS       (THD - BASED APPROACH)  ----- :  ", 10000 * (Cd_vrt_thd), "                                          |");
        printf("%s%8.1f%s\n",			"|                          ---------------------------------------------------- VORTEX FORCE  :  ", 10000 * (Cd_vf_thd), "     |");
        printf("%s%8.1f%s\n",			"|                          --------------------------------------- COMPRESSIBILITY CORRECTION :  ", 10000 * (Cd_compr), "     |");
        printf("%s%8.1f%s\n",			"|                          ------------------------------------- OUTER VORTICITY CONTRIBUTION :  ", 10000 * (Cd_fs_thd), "     |");
        printf("%s%8.1f%s\n",			"|                          --------------------------------------------- VISCOUS CONTRIBUTION :  ", 10000 * (Cd_fmu_thd), "     |");
        printf("%s%8.1f%s\n",			"|  TOTAL LIFT              (THD - BASED APPROACH)  ------ :  ", (Cl_vrt_thd), "                                         |");
        printf("%s%8.1f%s\n",			"|                          ---------------------------------------------------- VORTEX FORCE  :  ", (Cl_vf_thd), "     |");
        printf("%s%8.1f%s\n",			"|                          --------------------------------------- COMPRESSIBILITY CORRECTION :  ", (Cl_compr), "     |");
        printf("%s%8.1f%s\n",			"|                          ------------------------------------- OUTER VORTICITY CONTRIBUTION :  ", (Cl_fs_thd), "     |");
        printf("%s%8.1f%s\n",			"|                          --------------------------------------------- VISCOUS CONTRIBUTION :  ", (Cl_fmu_thd), "     |");
        printf("%s\n",					"|                                                                                                             |");
        printf("%s\n",					"|-------------------------------------------------------------------------------------------------------------|");  

        if (vrt_std_on == 1) {
            std::cout << "\n" << endl;
            printf("%s\n",					"|            FAR-FIELD AERODYNAMIC FORCE BREAK-DOWN  (LAMB VECTOR METHOD: WU ET AL. STANDARD APPROACH)        |");
            printf("%s\n",					"|------------------------------------------------------------------------------------------------------------ |");
            printf("%s%8.1f%s\n",	        "|                                                 X_w/c: ", (Xmax), "\t\t\t\t\t"                          "      |");
            printf("%s\n",					"|-------------------------------------------------------------------------------------------------------------|");
            printf("%s\n",					"|                                                                                                             |");
            printf("%s%8.1f%s\n",			"|  TOTAL DRAG COUNTS       (STD - APPROACH)  ----- :  ", 10000 * (Cd_vrt_std), "                                                |");
            printf("%s%8.1f%s\n",			"|                          ---------------------------------------------------- VORTEX FORCE  :  ", 10000 * (Cd_vf_std), "     |");
            printf("%s%8.1f%s\n",			"|                          --------------------------------------- COMPRESSIBILITY CORRECTION :  ", 10000 * (Cd_compr), "     |");
            printf("%s%8.1f%s\n",			"|                          ------------------------------------- OUTER VORTICITY CONTRIBUTION :  ", 10000 * (Cd_fs_std), "     |");
            printf("%s%8.1f%s\n",			"|                          --------------------------------------------- VISCOUS CONTRIBUTION :  ", 10000 * (Cd_fmu_std), "     |");
            printf("%s%8.1f%s\n",			"|  TOTAL LIFT              (STD - APPROACH)  ------ :  ", (Cl_vrt_std), "                                               |");
            printf("%s%8.1f%s\n",			"|                          ---------------------------------------------------- VORTEX FORCE  :  ", (Cl_vf_std), "     |");
            printf("%s%8.1f%s\n",			"|                          --------------------------------------- COMPRESSIBILITY CORRECTION :  ", (Cl_compr), "     |");
            printf("%s%8.1f%s\n",			"|                          ------------------------------------- OUTER VORTICITY CONTRIBUTION :  ", (Cl_fs_std), "     |");
            printf("%s%8.1f%s\n",			"|                          --------------------------------------------- VISCOUS CONTRIBUTION :  ", (Cl_fmu_std), "     |");
            printf("%s\n",					"|                                                                                                             |");
            printf("%s\n",					"|-------------------------------------------------------------------------------------------------------------|");  
        } 

        if (vrt_std_on == 1) {
            outputFile << Xmax << "\t\t" << 10000*Cd_PT << "\t\t" << 10000*Cd_wave_PT << "\t\t" << 10000*Cd_viscous_PT << "\t\t" << 10000*Cd_spurious_PT << "\t\t" << 10000*Cd_DV << "\t\t" << 10000*Cd_wave_DV << "\t\t" << 10000*Cd_viscous_DV << "\t\t" << 10000*Cd_spurious_DV << "\t\t"
            << 10000*Cd_vf_thd << "\t\t" << Cl_vf_thd << "\t\t" << 10000*Cd_vf_std << "\t\t" << Cl_vf_std << "\t\t" << 10000*Cd_compr << "\t\t" << Cl_compr << "\t\t" << 10000*Cd_fs_thd << "\t\t" << Cl_fs_thd << "\t\t" << 10000*Cd_fs_std << "\t\t" << Cl_fs_std << "\t\t" << 10000*Cd_fmu_thd << "\t\t" << Cl_fmu_thd << "\t\t" << 10000*Cd_fmu_std << "\t\t" << Cl_fmu_std << "\t\t" << 10000*Cd_vrt_thd << "\t\t" << Cl_vrt_thd << "\t\t" << 10000*Cd_vrt_std << "\t\t" << Cl_vrt_std << "\t\t" 
            << 10000*(F_s_x_thd / (q_inf*chord)) << "\t\t"  << 10000*(F_s_x_std / (q_inf*chord)) << "\t\t" << 10000*Cd_p_thd << "\t\t" << 10000*Cd_p_std << "\t\t" << 10000*Cd_i_thd << "\t\t" << 10000*Cd_i_std << "\t\t" << 10000*Cd_wv_thd << "\t\t" << 10000*Cd_vs_thd << "\t\t" << 10000*Cd_sp_thd << "\t\t" << 10000*Cd_wv_std << "\t\t" << 10000*Cd_vs_std << "\t\t" << 10000*Cd_sp_std << "\t\t" << 10000*Cd_nf << "\t\t" << Cl_nf << "\n";
        } else {
            outputFile << Xmax << "\t\t" << 10000*Cd_PT << "\t\t" << 10000*Cd_wave_PT << "\t\t" << 10000*Cd_viscous_PT << "\t\t" << 10000*Cd_spurious_PT << "\t\t" << 10000*Cd_DV << "\t\t" << 10000*Cd_wave_DV << "\t\t" << 10000*Cd_viscous_DV << "\t\t" << 10000*Cd_spurious_DV << "\t\t" << 10000*Cd_vf_thd << "\t\t" << Cl_vf_thd << "\t\t" << 10000*Cd_compr << "\t\t" << Cl_compr << "\t\t" << 10000*Cd_fs_thd << "\t\t" << Cl_fs_thd << "\t\t" << 10000*Cd_fmu_thd << "\t\t" << Cl_fmu_thd << "\t\t" << 10000*Cd_vrt_thd << "\t\t" << Cl_vrt_thd << "\t\t" 
            << 10000*(F_s_x_thd / (q_inf*chord)) << "\t\t" << "\t\t" << 10000*Cd_p_thd << "\t\t" << 10000*Cd_i_thd << "\t\t" << 10000*Cd_wv_thd << "\t\t" << 10000*Cd_sp_thd << "\t\t" << 10000*Cd_nf << "\t\t" << Cl_nf << "\n";
        }

        // Execution time for 1 iteration on Xw
        if(Xmax == 1){
            auto end_single = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_single = end_single - start_single;
            elapsed_time_single = elapsed_single.count();
        } 
    }

    // ************************** WRITING THE OUTPUT FILES (.szplt), (.txt) *****************************************

    // Making new variables recognaziable from the zone
    varTypes.resize(numVarsOut, FieldDataType_Double);
    shareVarFromZone.resize(numVarsOut, 0.0);
    passiveVarList.resize(numVarsOut, 0.0);
    valueLocation.resize(numVarsOut, 0.0); 

    // Add normals to cellCenteredVariables
    for(int cell = 0; cell < numCells; cell++) {
        for(int face = 0; face < 4; face++) {
            // Add normals to cellCenteredVariables in order to write these variables in output
            int nxID = 3 + face*2;
            int nyID = 4 + face*2;
            cellCenteredVariables[cell][nxID] = normals[cell][face*2]; //nx
            cellCenteredVariables[cell][nyID] = normals[cell][face*2 + 1]; //ny
        }
    }

    // New variables settings (ValueLocation = 0 -> node variable; ValueLocation = 1 -> cell centered variable)
    for (int32_t var = numVars; var < numVarsOut; var++) {
        varTypes[var] = FieldDataType_Double; // double type for all variables
    }

    for (int32_t var = 0; var < numVarsOut; var++) {
        passiveVarList[var] = 0; // no passive variabile
    }

    for (int32_t var = numVars; var < numVarsOut; var++) {
        if (var >= numVars && var < numVars + nodeVars) { // node variables
            valueLocation[var] = 1; 
        }
        else { // cell centered variables
            valueLocation[var] = 0;
        }
    } 

    for (int32_t var = 0; var < numVarsOut; var++) {
        if (var < numVars) {
            shareVarFromZone[var] = shareVarFromZone[var];
        }
        else {
            shareVarFromZone[var] = 0;
        }
    }
    
    // Opening the output file
        void* outputFileHandle = NULL;
        int32_t fileFormat = 1; // .szplt format
        R = tecFileWriterOpen(argv[3], dataSetTitle, varListStream.str().c_str(), fileFormat, 0, 1, NULL, &outputFileHandle);
        if (R != 0) throw std::runtime_error("Failed to open output file.");
        // R = tecFileSetDiagnosticsLevel(outputFileHandle, 1);
        tecStringFree(&dataSetTitle);
    
    // Creating the unstructured zone in the output file
        int32_t outputZone = 0;
        R = tecZoneCreateFE(outputFileHandle, zoneTitle, zoneType, iMax, numCells, &varTypes[0],
        &shareVarFromZone[0], &valueLocation[0], &passiveVarList[0],
        shareConnectivityFromZone, numFaceConnections, faceNeighborMode, &outputZone);

    // Writing nodes variables
     for (int32_t var = 1; var <= numVars + nodeVars; var++) {
         // Temporary vector to store nodes variables
         std::vector<double> tempValues(numNodes, 0.0);

         // Collect info from tempValues to nodeVariables
         for(size_t node = 0; node < numNodes; node++) {
             tempValues[node] = nodeVariables[node][var - 1];
             }

         // Writing
         R = tecZoneVarWriteDoubleValues(outputFileHandle, outputZone, var, 1, numNodes, &tempValues[0]);
        }

     // Writing connectivity informations
     if (zoneType != 0 && shareConnectivityFromZone == 0)
     {
         for (int64_t elem = 0; elem < numCells; elem++) {
         R = tecZoneNodeMapWrite32(outputFileHandle, outputZone, elem + 1, elem + 1, nodesPerCell, connectivity[elem].data());    
         }

     tecStringFree(&zoneTitle);
     }

     // Writing cell centered variables
    for (int32_t var = 1; var <= NUM_CELL_VARS; var++) {
        vector<double> tempValues(numCells, 0.0);
        for(size_t cell = 0; cell < numCells; cell++) {
            tempValues[cell] = cellCenteredVariables[cell][var - 1];
        }
        int32_t varID = numVars + nodeVars + var;
        R = tecZoneVarWriteDoubleValues(outputFileHandle, outputZone, varID, 1, numCells, &tempValues[0]);
    } 

    /*  // Writing cell centered TEMPORARY variables
    for (int32_t var = 1; var <= NUM_TEMP_CELL_VARS; var++) {
        vector<double> tempValues(numCells, 0.0);
        for(size_t cell = 0; cell < numCells; cell++) {
            tempValues[cell] = ccVariables[cell][var - 1];
        }
        int32_t varID = numVars + nodeVars + NUM_CELL_VARS + var;
        R = tecZoneVarWriteDoubleValues(outputFileHandle, outputZone, varID, 1, numCells, &tempValues[0]);
    }  */
    
    // Closing files
    tecFileReaderClose(&inputFileHandle);
    tecFileWriterClose(&outputFileHandle);
    outputFile.close();

    std::cout << "\n";
    auto end_total = std::chrono::high_resolution_clock::now(); // end
    std::chrono::duration<double> elapsed_total = end_total - start_total;
    std::cout << "PREP + EXECUTION TIME FOR ONE ITERATION ON Xw: " << elapsed_time_single << " s" << std::endl;
    std::cout << "TOTAL EXECUTION TIME: " << elapsed_total.count() << " s" << std::endl;
    std::cout << "EXECUTION COMPLETED :)" << std::endl;
    return 0;
}