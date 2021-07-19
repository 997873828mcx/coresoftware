// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4ZDCSTEPPINGACTION_H
#define G4DETECTORS_PHG4ZDCSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <Geant4/G4TouchableHandle.hh>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4ZDCDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class PHParameters;

class PHG4ZDCSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4ZDCSteppingAction(PHG4ZDCDetector*, const PHParameters* parameters);

  //! destroctor
  virtual ~PHG4ZDCSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  int FindIndex(G4TouchableHandle& touch, int& j, int& k);
  
  double ZDCResponce(double beta, double angle);
  
  double ZDCEResponce(double E, double angle);
 
  //! pointer to the detector
  PHG4ZDCDetector* m_Detector;

  //! pointer to hit container
  PHG4HitContainer* m_SignalHitContainer;
  PHG4HitContainer* m_AbsorberHitContainer;
  const PHParameters* m_Params;
  PHG4HitContainer* m_CurrentHitContainer;
  PHG4Hit* m_Hit;
  PHG4Shower* m_CurrentShower;

  int m_IsActiveFlag;
  int absorbertruth;
  int light_scint_model;
  int m_IsBlackHole;

  const std::array<std::array<double, 18>,9>m_PMMA05 = 
    {{
    {16.258, 6.771, 3.844, 2.432, 1.531, 1.003, 0.543, 0.195, 0.102,
     0.053, 0.017, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000} ,

    {14.722,10.130, 6.809, 4.436, 3.350, 2.558, 1.893, 1.404, 0.905,
     0.314, 0.173, 0.094, 0.046, 0.000, 0.000, 0.000, 0.000, 0.000}

    ,{9.922, 5.795, 6.179, 5.376, 4.152, 3.399, 2.779, 2.161, 1.745,
      1.171, 0.479, 0.233, 0.128, 0.080, 0.000, 0.000, 0.000, 0.000}

    ,{6.800, 3.826, 3.449, 4.751, 4.478, 3.885, 3.413, 2.713, 2.294,
      1.874, 1.264, 0.590, 0.301, 0.166, 0.095, 0.000, 0.000, 0.000}

    ,{4.222, 2.743, 2.408, 2.871, 4.146, 4.024, 3.620, 3.134, 2.713,
     2.232, 1.790, 1.187, 0.527, 0.284, 0.150, 0.076, 0.000, 0.000}

    ,{3.349, 2.428, 2.112, 2.526, 3.986, 4.061, 3.705, 3.200, 2.861,
      2.437, 1.790, 1.482, 0.589, 0.355, 0.196, 0.081, 0.000, 0.000}

    ,{1.825, 1.986, 1.961, 2.171, 3.075, 3.915, 3.783, 3.268, 2.912,
     2.561, 2.072, 1.659, 0.930, 0.425, 0.240, 0.131, 0.054, 0.000}

    ,{0.019, 1.612, 1.741, 1.827, 2.285, 3.659, 3.723, 3.512, 3.093,
     2.664, 2.319, 1.776, 1.283, 0.545, 0.328, 0.192, 0.084, 0.000}

    ,{0.000, 1.365, 1.446, 1.764, 2.082, 3.662, 3.548, 3.443, 3.128,
     2.777, 2.355, 1.936, 1.431, 0.636, 0.340, 0.235, 0.101, 0.000}

      } };
  const double m_Beta[9] = {0.70, 0.75, 0.80, 0.85, 0.90, 0.92, 0.95, 0.98, 1.00};
  const double m_BetaThersh = 0.671141;

  const std::array<std::array<double, 36>,11>m_PMMA05E =
     {{
	 {0.002, 0.009, 0.012, 0.047, 0.066, 0.053, 0.046, 0.041, 0.033,
	  0.014, 0.007, 0.005, 0.008, 0.001, 0.000, 0.000, 0.000, 0.000,
	  0.001, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
	 
	 {0.067, 0.141, 0.222, 0.319, 0.322, 0.296, 0.243, 0.212, 0.173,
	  0.113, 0.104, 0.052, 0.039, 0.025, 0.027, 0.008, 0.007, 0.007,
	  0.005, 0.004, 0.002, 0.005, 0.000, 0.001, 0.002, 0.000, 0.000,
	  0.001, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},

	 {0.213, 0.385, 0.607, 0.635, 0.652, 0.602, 0.529, 0.479, 0.386,
	  0.320, 0.235, 0.186, 0.135, 0.097, 0.074, 0.044, 0.031, 0.028,
	  0.020, 0.022, 0.006, 0.008, 0.015, 0.007, 0.003, 0.001, 0.004,
	  0.004, 0.000, 0.000, 0.002, 0.001, 0.002, 0.001, 0.000, 0.000},

	 {0.786, 1.278, 1.547, 1.718, 1.784, 1.748, 1.491, 1.442, 1.276,
	  1.088, 0.839, 0.726, 0.473, 0.385, 0.275, 0.200, 0.165, 0.103,
	  0.079, 0.074, 0.069, 0.049, 0.054, 0.021, 0.026, 0.024, 0.016,
	  0.016, 0.009, 0.014, 0.011, 0.019, 0.008, 0.008, 0.001, 0.000},

	 {1.622, 2.173, 2.584, 2.687, 2.892, 2.845, 2.709, 2.448, 2.254,
	  1.913, 1.579, 1.201, 0.883, 0.609, 0.390, 0.252, 0.147, 0.133,
	  0.071, 0.067, 0.042, 0.023, 0.029, 0.022, 0.006, 0.018, 0.012,
	  0.010, 0.004, 0.009, 0.006, 0.006, 0.004, 0.005, 0.007, 0.000},

	 {1.984, 2.660, 3.081, 3.080, 3.245, 3.449, 3.299, 3.096, 2.906,
	  2.411, 2.063, 1.574, 1.049, 0.680, 0.411, 0.236, 0.097, 0.047,
	  0.025, 0.017, 0.007, 0.006, 0.002, 0.006, 0.001, 0.008, 0.001,
	  0.007, 0.005, 0.001, 0.001, 0.005, 0.015, 0.006, 0.009, 0.000},

	 {1.884, 2.399, 2.716, 2.779, 3.147, 3.487, 3.497, 3.365, 3.043,
	  2.694, 2.232, 1.730, 1.186, 0.684, 0.359, 0.207, 0.073, 0.021,
	  0.017, 0.017, 0.011, 0.000, 0.002, 0.001, 0.001, 0.002, 0.011,
	  0.000, 0.003, 0.000, 0.008, 0.006, 0.014, 0.002, 0.000, 0.000},
	
	 {1.567, 1.873, 2.094, 2.311, 2.704, 3.493, 3.590, 3.454, 3.271,
	  2.716, 2.385, 1.909, 1.421, 0.591, 0.353, 0.209, 0.112, 0.007,
	  0.004, 0.000, 0.004, 0.001, 0.000, 0.004, 0.002, 0.004, 0.000,
	  0.000, 0.000, 0.003, 0.000, 0.000, 0.005, 0.004, 0.000, 0.000},
	 
	 {1.142, 1.426, 1.583, 1.928, 2.336, 3.395, 3.594, 3.382, 3.180,
	  2.907, 2.347, 1.961, 1.445, 0.613, 0.382, 0.232, 0.109, 0.003,
	  0.000, 0.000, 0.000, 0.000, 0.000, 0.001, 0.000, 0.001, 0.000,
	  0.001, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},

	 {0.638, 1.437, 1.560, 1.646, 2.106, 3.539, 3.864, 3.497, 3.223,
	  2.732, 2.345, 1.960, 1.440, 0.680, 0.367, 0.188, 0.099, 0.002,
	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},

	 {0.020, 1.407, 1.434, 1.720, 2.081, 3.498, 3.577, 3.495, 3.205,
	  2.788, 2.350, 2.039, 1.428, 0.684, 0.372, 0.201, 0.119, 0.000,
	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000}
       }}; 
   const double m_E[11] = {0.0005,0.00055,0.0006,0.00075,0.001,0.0015, 0.002,0.003,0.005,0.01,0.05};
  
};

#endif  // G4DETECTORS_PHG4ZDCSTEPPINGACTION_H
