# single crystal one element calibration
# for the HR-DIC experiment on 316L stainless steel

[GlobalParams]
  displacements = 'ux uy uz'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmax = 1.0
  ymax = 1.0
  zmax = 1.0
  elem_type = HEX8
  displacements = 'ux uy uz'
[]

[Variables]
  [./ux]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]

  [./uy]
    order = FIRST
    family = LAGRANGE
  [../]

  [./uz]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Modules/TensorMechanics/Master/all]
  strain = FINITE
  add_variables = true
  generate_output = stress_xx
[]

[AuxVariables]

  # slip_increment is the slip rate
  # so units are 1/time
  [./slip_increment_1]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_2]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_3]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_4]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_5]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_6]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_7]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_8]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_9]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_10]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_11]
   order = FIRST
   family = MONOMIAL
  [../]
  [./slip_increment_12]
   order = FIRST
   family = MONOMIAL
  [../]

  [./slip_increment_vector]
    order = FIRST
    family = MONOMIAL
    components = 12
  [../]

  [./dslip_increment_dedge]
    order = CONSTANT
    family = MONOMIAL
    components = 12
  [../]
  
  [./dslip_increment_dscrew]
    order = CONSTANT
    family = MONOMIAL
    components = 12
  [../]
  
  [./rho_ssd_1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./slip_resistance_1]
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[UserObjects]
  [./prop_read]
    type = GrainPropertyReadFile
    prop_file_name = 'EulerOneElement.txt'
    nprop = 3
    ngrain = 1
    read_type = indexgrain
  [../]
[]

[Functions]

  [./disp_load]
    type = ParsedFunction
    value = '0.2230*t/60.0'
  [../]

[]

[AuxKernels]

  [./slip_increment_1]
   type = MaterialStdVectorAux
   variable = slip_increment_1
   property = slip_increment
   index = 0
   execute_on = timestep_end
  [../]
  [./slip_increment_2]
   type = MaterialStdVectorAux
   variable = slip_increment_2
   property = slip_increment
   index = 1
   execute_on = timestep_end
  [../]
  [./slip_increment_3]
   type = MaterialStdVectorAux
   variable = slip_increment_3   
   property = slip_increment
   index = 2
   execute_on = timestep_end
  [../]
  [./slip_increment_4]
   type = MaterialStdVectorAux
   variable = slip_increment_4
   property = slip_increment
   index = 3
   execute_on = timestep_end
  [../]
  [./slip_increment_5]
   type = MaterialStdVectorAux
   variable = slip_increment_5
   property = slip_increment
   index = 4
   execute_on = timestep_end
  [../]
  [./slip_increment_6]
   type = MaterialStdVectorAux
   variable = slip_increment_6
   property = slip_increment
   index = 5
   execute_on = timestep_end
  [../]
  [./slip_increment_7]
   type = MaterialStdVectorAux
   variable = slip_increment_7   
   property = slip_increment
   index = 6
   execute_on = timestep_end
  [../]
  [./slip_increment_8]
   type = MaterialStdVectorAux
   variable = slip_increment_8
   property = slip_increment
   index = 7
   execute_on = timestep_end
  [../]
  [./slip_increment_9]
   type = MaterialStdVectorAux
   variable = slip_increment_9
   property = slip_increment
   index = 8
   execute_on = timestep_end
  [../]
  [./slip_increment_10]
   type = MaterialStdVectorAux
   variable = slip_increment_10
   property = slip_increment
   index = 9
   execute_on = timestep_end
  [../]
  [./slip_increment_11]
   type = MaterialStdVectorAux
   variable = slip_increment_11   
   property = slip_increment
   index = 10
   execute_on = timestep_end
  [../]
  [./slip_increment_12]
   type = MaterialStdVectorAux
   variable = slip_increment_12
   property = slip_increment
   index = 11
   execute_on = timestep_end
  [../]

  [./build_slip_increment_vector]
    type = BuildArrayVariableAux
    variable = slip_increment_vector
    component_variables = 'slip_increment_1 slip_increment_2 slip_increment_3 slip_increment_4 slip_increment_5 slip_increment_6 slip_increment_7 slip_increment_8 slip_increment_9 slip_increment_10 slip_increment_11 slip_increment_12'
  [../]

  [./edge_directional_derivative]
    type = ArrayDirectionalDerivative
    variable = dslip_increment_dedge
    gradient_variable = slip_increment_vector
    dislo_character = edge
  	execute_on = timestep_end
  [../]
  
  [./screw_directional_derivative]
    type = ArrayDirectionalDerivative
    variable = dslip_increment_dscrew
    gradient_variable = slip_increment_vector
    dislo_character = screw
  	execute_on = timestep_end
  [../]  
  
  [./rho_ssd_1]
    type = MaterialStdVectorAux
    variable = rho_ssd_1
    property = rho_ssd
    index = 0
    execute_on = timestep_end
  [../]
  
  [./slip_resistance_1]
    type = MaterialStdVectorAux
    variable = slip_resistance_1
    property = slip_resistance
    index = 0
    execute_on = timestep_end
  [../]
  
[]

[BCs]
  [./z0_back]
    type = DirichletBC
    variable = uz
    boundary = back
    value = 0.0
  [../]

  [./y0_bottom]
    type = DirichletBC
    variable = uy
    boundary = bottom
    value = 0.0
  [../]
  
  [./x0_left]
    type = DirichletBC
    variable = ux
    boundary = left
    value = 0.0
  [../]

  [./x1_right]
    type = FunctionDirichletBC
    variable = ux
    boundary = right
    function = disp_load
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCPGrain
    C_ijkl = '2.046e5 1.377e5 1.377e5 2.046e5 1.377e5 2.046e5 1.262e5 1.262e5 1.262e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
  [../]
  [./stress]
    type = ComputeDislocationCrystalPlasticityStress
    crystal_plasticity_models = 'trial_xtalpl'
    tan_mod_type = exact
    maximum_substep_iteration = 1
  [../]
  [./trial_xtalpl]
    type = CrystalPlasticityDislocationUpdate
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.txt
	ao = 0.001
	xm = 0.1
	burgers_vector_mag = 0.000256
	shear_modulus = 86000.0 # MPa
	alpha_0 = 0.3
	r = 1.4
	tau_c_0 = 0.112
	k_0 = 0.053000000000000005
	y_c = 0.0026
	init_rho_ssd = 20.5
	init_rho_gnd_edge = 0.0
	init_rho_gnd_screw = 0.0
	# These activate slip gradients
	# they are compulsory
	# codes currently has problems if not introduced
	# to remove slip gradients, zero arrays can be passed
	dslip_increment_dedge = dslip_increment_dedge
	dslip_increment_dscrew = dslip_increment_dscrew
  [../]
[]

[Postprocessors]
  [stress_on_load_surface]
    type = ElementAverageValue
    variable = stress_xx
	block = 0
    use_displaced_mesh = true
  []
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.05
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -ksp_type -ksp_gmres_restart'
  petsc_options_value = ' asm      2              lu            gmres     200'
  nl_abs_tol = 1e-8
  nl_rel_step_tol = 1e-8
  dtmax = 0.05
  nl_rel_tol = 1e-8
  end_time = 60.0
  dtmin = 0.00001
  #num_steps = 100 #100 to see plastic shear stress
  nl_abs_step_tol = 1e-8
[]

[Outputs]
  exodus = true
  csv = true
  interval = 10
[]
