#ifndef NEWLOGLIKFITTER_PREMPSFITDRIVER_H
#define NEWLOGLIKFITTER_PREMPSFITDRIVER_H


void newLogLikFitter_preMPSfitdriver(
    const std::string &fit_description_string,
    const std::string &oname_prepend,
    const std::string &oname_append,
    const std::string &odir,
    double *set_min_point, double &set_min_point_fval)
{

    // initialize
    std::cout << "seo=" << gSystematics.systematic_energy_offset << std::endl;
    std::cout << "sem=" << gSystematics.systematic_energy_scale << std::endl;
    std::cout << "sef=" << gSystematics.systematic_efficiency << std::endl;
    std::cout << "sen=" << gSystematics.systematic_enrichment << std::endl;
    std::cout << "seos=" << gSystematics.systematic_energy_offsetsmall << std::endl;
    std::cout << "sftV=" << gSystematics.systematic_foil_thickness_virtual << std::endl;
    std::cout << "selV=" << gSystematics.systematic_dEdX_virtual << std::endl;
    std::cout << "sbV=" << gSystematics.systematic_brem_virtual << std::endl;
    std::cout << "sftN=" << gSystematics.systematic_foil_thickness_nominal << std::endl;
    std::cout << "selN=" << gSystematics.systematic_dEdX_nominal << std::endl;
    std::cout << "sbN=" << gSystematics.systematic_brem_nominal << std::endl;
    std::string name_extra_g_mode_fake_data;
    if(g_mode_fake_data == true)
    {
        name_extra_g_mode_fake_data = "fakedata";
    }
    else
    {
        name_extra_g_mode_fake_data = "data";
    }
    std::string name_extra = name_extra_g_mode_fake_data
                           + "_seo_" + std::to_string(gSystematics.systematic_energy_offset)
                           + "_sem_" + std::to_string(gSystematics.systematic_energy_scale)
                           + "_sef_" + std::to_string(gSystematics.systematic_efficiency)
                           + "_sen_" + std::to_string(gSystematics.systematic_enrichment)
                           + "_seos_" + std::to_string(gSystematics.systematic_energy_offsetsmall)
                           + "_sftV_" + std::to_string(gSystematics.systematic_foil_thickness_virtual)
                           + "_selV_" + std::to_string(gSystematics.systematic_dEdX_virtual)
                           + "_sbV_" + std::to_string(gSystematics.systematic_brem_virtual)
                           + "_sftN_" + std::to_string(gSystematics.systematic_foil_thickness_nominal)
                           + "_selN_" + std::to_string(gSystematics.systematic_dEdX_nominal)
                           + "_sbN_" + std::to_string(gSystematics.systematic_brem_nominal);

    // assuming that xi31 and 150Nd amplitude are free
    // this may break if the parameter_names.lst file is changed
    gNumberFreeParams = 2;

    // create minimizer
    ROOT::Minuit2::MnUserParameterState theParameterStateBefore;
    ROOT::Minuit2::VariableMetricMinimizer theMinimizer;
    MinimizeFCNAxialVector theFCN;

    // initialize fit
    //fitBackgrounds_init(theParameterState, theMinimizer, AdjustActs, AdjustActs_Err);
    const int xi_31_param_number = g_pg.get_xi_31_ext_param_number();
    const double xi_31_value = g_pg.file_params.at(xi_31_param_number).paramInitValue;
    const double xi_31_error = g_pg.file_params.at(xi_31_param_number).paramInitError;
    std::cout << "xi_31_param_number=" << xi_31_param_number
              << " xi_31=" << xi_31_value << " +- " << xi_31_error << std::endl;
    fitBackgrounds_init(theParameterStateBefore, theMinimizer, xi_31_value, xi_31_error);

    // get parameters and chi2 value before fit
    std::vector<double> params_before = theParameterStateBefore.Params();
    std::vector<double> param_errs_before = theParameterStateBefore.Errors();
    double fval_before = theFCN.operator()(params_before);
    //int ndf = theFCN.ndf - theParameterStateBefore.VariableParameters();
    //int ndf = theFCN.nch - g_pg.get_number_free_params();
    int nch = theFCN.nch;
    //int nfp = g_pg.get_number_free_params();
    int nfp = gNumberFreeParams;
    int ndf = nch - nfp;

    std::cout << "nch=" << theFCN.nch << std::endl;
    std::cout << "gNumberFreeParams=" << gNumberFreeParams << std::endl;
    //std::cout << "g_pg.get_number_free_params()=" << g_pg.get_number_free_params() << std::endl;

/*
    std::cout << "fval=" << fval_before << std::endl;
    for(int i = 0; i < params_before.size(); ++ i)
    {
        std::cout << "i=" << i << " params_before[" << i << "]=" << params_before[i] << std::endl;
    }
*/

    // draw before fit
    draw_input_data drawinputdata;
    drawinputdata.chi2 = fval_before;
    drawinputdata.nch = nch;
    drawinputdata.nfp = nfp;
    //drawinputdata.serial_dir = "xifree";
    drawinputdata.serial_dir = odir;
    //drawinputdata.saveas_filename = std::string("xifree_before") + "_" + name_extra;
    drawinputdata.saveas_filename = oname_prepend + "before" + oname_append;
    drawinputdata.saveas_png = true;
   
    draw(drawinputdata,
         params_before,
         param_errs_before);

    // exec fit
    // do fit with all parameters free
    ROOT::Minuit2::FunctionMinimum FCN_min =
        fitBackgrounds_exec(
            theParameterStateBefore,
            theMinimizer,
            theFCN);

    // get result
    ROOT::Minuit2::MnUserParameterState theParameterStateAfter = FCN_min.UserParameters();
    std::vector<double> params_after = theParameterStateAfter.Params();
    std::vector<double> param_errs_after = theParameterStateAfter.Errors();

    double fval_after = theFCN.operator()(params_after);
    //ndf = theFCN.ndf - theParameterStateAfter.VariableParameters();
    //ndf = theFCN.ndf - g_pg.get_number_free_params();
    nch = theFCN.nch;
    nfp = gNumberFreeParams;
    ndf = nch - nfp;

/*
    std::cout << "fval=" << fval_after << std::endl;
    for(int i = 0; i < params_after.size(); ++ i)
    {
        std::cout << "i=" << i << " params_after[" << i << "]=" << params_after[i] << std::endl;
    }
*/

    // draw result
    drawinputdata.chi2 = fval_after;
    drawinputdata.nfp = nfp;
        //drawinputdata.nch = nch; // this could probably change in theory
        //drawinputdata.nfp = nfp; // these probably do not change
    drawinputdata.serial_dir = odir;
    //drawinputdata.saveas_filename = std::string("xifree_after") + "_" + name_extra;
    drawinputdata.saveas_filename = oname_prepend + "after" + oname_append;
   
    draw(drawinputdata,
         params_after,
         param_errs_after);

    
    // minimize
    //ROOT::Minuit2::FunctionMinimum FCN_min = theMinimizer.Minimize(theFCN, init_par, init_err);
    //ROOT::Minuit2::FunctionMinimum FCN_min = theMinimizer.Minimize(theFCN, init_par, init_err);
    /*
    std::cout << "Minimization finished" << std::endl;
    std::cout << "minimum: " << FCN_min << std::endl;
    std::cout << "chi2: " << FCN_min.Fval() << std::endl;
    std::cout << "edm: " << FCN_min.Edm() << std::endl;
    */


    //std::cout << "fval_after=" << fval_after << " for params_after[0]=" << params_after[0] << " params_after[1]=" << params_after[1] << std::endl;
    //std::cout << "fval_before=" << fval_before << std::endl;

    set_min_point[0] = params_after.at(1);
    set_min_point[1] = params_after.at(0);
    set_min_point_fval = fval_after;

    std::cout << fit_description_string << std::endl;
    std::cout << "SYSTEMATICS: CONSTANT OFFSET: " << gSystematics.systematic_energy_offset << " MeV" << std::endl;
    std::cout << "SYSTEMATICS: CONSTANT SCALE: " << gSystematics.systematic_energy_scale << " MeV" << std::endl;
    std::cout << "SYSTEMATICS: EFFICIENCY: " << gSystematics.systematic_efficiency << "" << std::endl;
    std::cout << "SYSTEMATICS: ENRICHMENT: " << gSystematics.systematic_enrichment << "" << std::endl;
    std::cout << "SYSTEMATICS: OFFSETSMALL: " << gSystematics.systematic_energy_offsetsmall << " MeV" << std::endl;
    std::cout << "SYSTEMATICS: FOIL THICKNESS (V): " << gSystematics.systematic_foil_thickness_virtual << " " << std::endl;
    std::cout << "SYSTEMATICS: ENERGY LOSS (V): " << gSystematics.systematic_dEdX_virtual << " " << std::endl;
    std::cout << "SYSTEMATICS: BREMSSTRAHLUNG (V): " << gSystematics.systematic_brem_virtual << " " << std::endl;
    std::cout << "SYSTEMATICS: FOIL THICKNESS (N): " << gSystematics.systematic_foil_thickness_nominal << " " << std::endl;
    std::cout << "SYSTEMATICS: ENERGY LOSS (N): " << gSystematics.systematic_dEdX_nominal << " " << std::endl;
    std::cout << "SYSTEMATICS: BREMSSTRAHLUNG (N): " << gSystematics.systematic_brem_nominal << " " << std::endl;
    std::cout << "Result: " << std::endl;
    std::cout << "fval_before=" << fval_before << std::endl;
    std::cout << "fval_after=" << fval_after
              << " for params_after[0]=" << params_after[0] << " +- " << param_errs_after[0]
              << " params_after[1]=" << params_after[1] << " +- " << param_errs_after[1] << std::endl;
    std::cout << std::endl;
}


#endif // NEWLOGLIKFITTER_PREMPSFITDRIVER_H
