function run_core_tests()
%RUN_CORE_TESTS Run self-contained public regressions (no measured data needed).
verify_frft;
verify_autofocus;
verify_geometry_output;
verify_adaptive_constraint;
verify_project_entrypoints;
verify_workflow_profiles;
fprintf('All core tests passed.\n');
end
