function run_all_tests()
%RUN_ALL_TESTS Run the deterministic smoke suite in implementation order.

run_structure_topology_tests();
run_static_tests();
run_rhs_tests();
run_integrator_tests();
run_ground_state_tests();
run_interface_tests();
run_open_sot_tests();
fprintf('All CrI3 smoke tests passed.\n');
end
