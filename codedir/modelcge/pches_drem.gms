$title PCHES Dynamic Regional Economic Model - Open Computable General Equilibrium Model with Coupled Integration, Intra-national Trade, and Differentiated Electricity Production

$ontext
This dynamic CGE model is designed for coupling with the power system model. We incorperate
the Sullivan productivity damage function is integrated at each time period, but in the current version,
it is only active in the final period.
$offtext

* ----------------------------------------------------------------------------*
*       Set system options:
* ----------------------------------------------------------------------------*

*               Decrease the size of the listing file
$offsymxref
$offsymlist

*               MPSGE Output options
*option sysout = on;

* ---------------------------------------------------------------------------- *
*       Select data target names:
* ---------------------------------------------------------------------------- *

*               Main Model Data Target Parameters

* Number of regions
$if not set t_r $set t_r 12

* Number of sectors
$if not set t_s $set t_s 30

* Last year
$if not set t_y $set t_y 2050
* SSP
$if not set t_p $set t_p 2

* GCM for Demand Response
$if not set t_g $set t_g GFDL-CM3

* Switch targets for time period and tfp calibration
$if not set s_tfp $set s_tfp 0

* ---------------------------------------------------------------------------- *
*       Select coupled model:
* ---------------------------------------------------------------------------- *

*               Switch targets for coupling

* WDRG
$if not set s_wd $set s_wd 0

* PSM (no productivity impacts)
$if not set s_pw $set s_pw 1

* PSM (with productivity impacts)
$if not set s_pf $set s_pf 1

* Demand response
$if not set s_dr $set s_dr 0

* ---------------------------------------------------------------------------- *
*       Define data based on target names:
* ---------------------------------------------------------------------------- *

*               Data targets

* IMPLAN Data
$if not set target $set target implan440-%t_s%sector-%t_r%region

* Population and population growth rate
$set target_pop popproj-SSP%t_p%-2100
$set target_pop_rate popproj-cagr-SSP%t_p%-2100

* ---------------------------------------------------------------------------- *
*       Define export name modification name:
* ---------------------------------------------------------------------------- *

*               Export targets
$if not set mod $set mod test_manual

* ---------------------------------------------------------------------------- *
*       Select coupled model data targets:
* ---------------------------------------------------------------------------- *

* WDRG
$if not set wdrg_run $set wdrg_run state_scenario1_iteration2

* Power System
$if not set psm_sul $set psmsul MtoCGE
$if not set psm_tfp $set psmtfp MtoCGE2

* Demand response does not need specifice output files.

* ----------------------------------------------------------------------------*
*       Select OS:
* ----------------------------------------------------------------------------*

$if not set s_os $set s_os 1

* ----------------------------------------------------------------------------*
*       Create folders for save results:
* ----------------------------------------------------------------------------*

*               Create under Windows

$if not %s_os%==0 $goto folder_unix

* Create results folders if necessary
$call 'if not exist output\nul mkdir output\'
$call 'if not exist output\%target%\nul mkdir output\%target%\'
$call 'if not exist listings\nul mkdir listings\'
$call 'if not exist logs\nul mkdir logs\'

* Create main results directory
$call 'if not exist output\%target%\gdx\nul mkdir output\%target%\gdx\'
$call 'if not exist output\%target%\csv\nul mkdir output\%target%\csv\'

* Create coupling directory for output
$call 'if not exist output\%target%\coupling\nul mkdir output\%target%\coupling\'

* Create WDRG coupling directory
$call 'if not exist output\%target%\coupling\wdrg\nul mkdir output\%target%\coupling\wdrg\'

* Create demand response coupling directory
$call 'if not exist output\%target%\coupling\drm\nul mkdir output\%target%\coupling\drm\'

* Create PSM coupling directory
$call 'if not exist output\%target%\coupling\psm\nul mkdir output\%target%\coupling\psm\'


*               Create under Unix

$label folder_unix
$if not %s_os%==0 $goto load_win

* Create results folders if necessary
$call 'if not exist output/nul mkdir output/'
$call 'if not exist listings/nul mkdir listings/'
$call 'if not exist logs/nul mkdir logs/'

* Create main results directory
$call 'if not exist output/gdx/nul mkdir output/gdx/'
$call 'if not exist output/csv/nul mkdir output/csv/'

* Create coupling directory for output
$call 'if not exist output/coupling/nul mkdir output/coupling/'

* Create WDRG coupling directory
$call 'if not exist output/coupling/wdrg/nul mkdir output/coupling/wdrg/'

* Create demand response coupling directory
$call 'if not exist output/coupling/drm/nul mkdir output/coupling/drm/'

* Create PSM coupling directory
$call 'if not exist output/coupling/psm/nul mkdir output/coupling/psm/'

* ---------------------------------------------------------------------------- *
*       Define other dataset locations:
* ---------------------------------------------------------------------------- *

*               Create under Windows

$label load_win
$if not %s_os%==0 $goto load_unix

* Load IMPLAN
$set load_set      "models\regionaldata.gms"
$set load_implan   'data\%target%\%target%_dtrdbal.gdx'
$set load_map      "data\%target%\agg_mapping.gms"
$set load_map2      "data\%target%\agg_mapping_psm.gms"

$label implan_wdrg_data_win
$if not %s_wd%==1 $goto next_load_win
*$set load_set      "models\regionaldata.gms"
*$set load_implan   'data\%target%\%target%_dtrdbal.gdx'
*$set load_map      "data\%target%\agg_mapping.gms"
$set load_wdrg_pre "models\process_inputs.gms"
$set load_yield_in 'data\wdrg\yield_%wdrg_run%.dat'
$set load_gcam     "data\wdrg\gcam_land_%wdrg_run%.dat"
$set load_yield_ot 'data\wdrg\yield_change_%wdrg_run%'

* Load elasticities
$label next_load_win
$set load_es       'data\%target%\es_drem.gdx'

* Load population
$set load_pop      'data\%target%\population\%target_pop%.gdx'
$set load_pop_carg 'data\%target%\population\%target_pop_rate%.gdx'

* Load macro parameters
$set load_macro    'data\%target%\macro_parameters.gdx'

* Load demand response data
$set load_drm      'data\%target%\demand-response-emulator\drm-%t_g%.gdx'

* Load psm electicity data
$set load_psm_ele  'data\power-system-model\%psmtfp%'

* Load psm productivity data
$set load_psm_nele 'data\power-system-model\%psmsul%'

*               Create under Unix

$label load_unix
$if not %s_os%==1 $goto results_win

* Load IMPLAN
$set load_set      "regionaldata.gms"
$set load_implan   'data/%target%/%target%_dtrdbal.gdx'
$set load_map      "data/%target%/agg_mapping.gms"
$set load_map2      "data/%target%/agg_mapping_psm.gms"

$label wdrg_load_unix
$if not %s_wd%==1 $goto next_load_unix
*$set load_set      "regionaldata.gms"
*$set load_implan   'data/%target/%target%_dtrdbal.gdx'
*$set load_map      "data/%target%/agg_mapping.gms"
$set load_wdrg_pre "process_inputs.gms"
$set load_yield_in 'data/wdrg/yield_%wdrg_run%.dat'
$set load_gcam     "data/wdrg/gcam_land_%wdrg_run%.dat"
$set load_yield_ot 'data/wdrg/yield_change_%wdrg_run%'

* Load elasticities
$label next_load_unix
$set load_es       'data/%target%/es_drem.gdx'

* Load population
$set load_pop      'data/%target%/population/%target_pop%.gdx'
$set load_pop_carg 'data/%target%/population/%target_pop_rate%.gdx'

* Load macro parameters
$set load_macro    'data/%target%/macro_parameters.gdx'

* Load demand response data
$set load_drm      'data/%target%/demand-response-emulator/drm-%t_g%.gdx'

* Load psm electicity data
$set load_psm_ele  '%psmtfp%'

* Load psm productivity data
$set load_psm_nele '%psmsul%'

* ---------------------------------------------------------------------------- *
*       Define result save locations:
* ---------------------------------------------------------------------------- *

*               Save results under Windows

$label results_win
$if not %s_os%==0 $goto results_unix

* Save level values
$set save_valuelevel 'output\%target%\gdx\%mod%_valuelevel.gdx'

* Save percentage change values
$set save_valuepchg  'output\%target%\gdx\%mod%_valuepchg.gdx'

* Save growth rates
$set save_growthrate 'output\%target%\gdx\%mod%_growthrate.gdx'

* Save model parameters
$set save_parameters 'output\%target%\gdx\%mod%_parameters.csv'

$set save1_parameters 'output\%target%\csv\%mod%_valuelevel.csv'
$set save2_parameters 'output\%target%\csv\%mod%_valuepchg.csv'
$set save3_parameters 'output\%target%\csv\%mod%_growthrate.csv'
*$set save4_parameters 'output\%target%\csv\%mod%_'
*$set save5_parameters 'output\%target%\csv\%mod%_'
*$set save6_parameters 'output\%target%\csv\%mod%_'

* Save WDRG coupling
$set save_wdrg_couple 'output\%target%\coupling\wdrg\DREM_output_%wdrg_run%.gdx'
$set save_wdrg_all    'output\%target%\coupling\wdrg\DREM_%wdrg_run%.gdx'

* Save demand coupling
$set save_drm        'output\%target%\coupling\drm\drm_output.gdx'

* Save psm coupling
$set save_psm_couple 'output\%target%\coupling\psm\CGEtoM.gdx'
$set save_psm_all    'output\%target%\coupling\psm\CGEtoM2.gdx'


*               Save results under Unix

$label results_unix
$if not %s_os%==1 $goto target_data

* Save level values
$set save_valuelevel 'output/gdx/%mod%_level'

* Save percentage change values
$set save_valuepchg  'output/gdx/%mod%_pchg'

* Save growth rates
$set save_growthrate 'output/gdx/%mod%_growth'

* Save model parameters
$set save_parameters 'output/gdx/%mod%_parameters.csv'

$set save1_parameters 'output/csv/%mod%_level.csv'
$set save2_parameters 'output/csv/%mod%_pchg.csv'
$set save3_parameters 'output/csv/%mod%_growth.csv'
*$set save4_parameters 'output/csv/%mod%_'
*$set save5_parameters 'output/csv/%mod%_'
*$set save6_parameters 'output/csv/%mod%_'

* Save WDRG coupling
$set save_wdrg_couple 'output/coupling/wdrg/REM_output_%wdrg_run%.gdx'
$set save_wdrg_all    'output/coupling/wdrg/DREM_%wdrg_run%.gdx'

* Save demand coupling
$set save_drm        'output/coupling/drm/%mod%_drm_output'

* Save psm coupling
$set save_psm_couple 'CGEtoM'
$set save_psm_all    'CGEtoM2'

* ---------------------------------------------------------------------------- *
*       Read the IMPLAN target dataset:
* ---------------------------------------------------------------------------- *

$label target_data

$include %load_set%

* ---------------------------------------------------------------------------- *
*       Set flags for calling coupled models:
* ---------------------------------------------------------------------------- *


PARAMETER
         switch_time     "Switch to select 2010-2050"
         switch_wdrg     "Switch to turn on WDRG coupling: == 1 on; == 0 off"
         switch_psm      "Switch to turn on PSM electricity impacts coupling: == 1 on; == 0 off"
         switch_psmf     "Switch to turn on PSM w/Sullivan impacts coupling: == 1 on; == 0 off"
         switch_dre      "Switch to turn on demand response emulator coupling: == 1 on; == 0 off"
;

switch_time = %t_y%;
switch_wdrg = %s_wd%;
switch_psm  = %s_pw%;
switch_psmf = %s_pf%;
switch_dre  = %s_dr%;

* ---------------------------------------------------------------------------- *
*       Set flags for calling tfp calibration values:
* ---------------------------------------------------------------------------- *

PARAMETER
         switch_tfp     "Switch to turn on tfp calibration: == 1 on; == 0 off"
;

switch_tfp = %s_tfp%;

* ---------------------------------------------------------------------------- *
*       Install parameters for couple model and TFP:
* ---------------------------------------------------------------------------- *

PARAMETER
        tfp(r,s)                Total factor productivity - region r for sector s
        wdrgtfp(r,s)            WDRG Total factor productivity - region r for sector s
        dre_prod(r,g,s)         Demand response emulator productivity - region r for good electricity and sector s
        dre_hh(r,s,h)           Demand response emulator productivity - region r for good electricity and household h
        psmtfp(r,g,s)           PSM productivity shock for sector s (electricity and Sullivan)
        psmtfp_fa(r,s)          PSM productivity shock for factors fa (electricity and Sullivan)
;

tfp(r,s)        = 1.00;
wdrgtfp(r,s)    = 1.00;
dre_prod(r,g,s) = 1.00;
dre_hh(r,s,h)   = 1.00;
psmtfp(r,g,s)   = 1.00;
psmtfp_fa(r,s)  = 1.00;

* ---------------------------------------------------------------------------- *
*       Define factor, good, and sector subsets:
* ---------------------------------------------------------------------------- *

SET
* Factor subsets
        fa(f)           Capital and labor factors          / empl, prop, othp /
        ft(f)           Business tax factor                / btax  /
* Sector subsets
        agrs(s)         Agricultural sectors               / apa, cba, frs, grn, oca,
                                                             osa,  pfb, vna /
        crps(s)         Crop sectors                       / cba, grn, oca, osa, pfb,
                                                             vna /
        rexs(s)         Resource extraction sector         / coa, ngd, cru, min /
        mats(s)         Materials sector                   / con, ele, per, trn
                                                             bom, cem, cpm, fbm,
                                                             pfm, tec, trm, wpm,
                                                             bos, fin, hlt, pub,
                                                             rtl, tel  /
* Goods subsets for production
        agrg(g)         Agricultural goods                 / apa, cba, frs, grn, oca,
                                                             osa,  pfb, vna /
        matg(g)         Material goods                     / apa, cba, frs, grn, oca,
                                                             osa,  pfb, vna,
                                                             con,
                                                             bom, cem, cpm, fbm, pfm,
                                                             tec, trm, wpm,
                                                             bos, fin, hlt, pub, rtl,
                                                             tel,
                                                             trn /
        eleg(g)         Electricity goods                  / ele /
        peng(g)         Primary energy goods               / coa, cru, ngd, per /
* secort subsets for consumption
        agrc(s)         Agricultural goods                 / apa, cba, frs, grn, oca,
                                                             osa,  pfb, vna /
        matc(s)         Material goods                     / apa, cba, frs, grn, oca,
                                                             osa,  pfb, vna,
                                                             con,
                                                             bom, cem, cpm, fbm, pfm,
                                                             tec, trm, wpm,
                                                             bos, fin, hlt, pub, rtl,
                                                             tel,
                                                             trn /
        elec(s)         Electricity goods                  / ele /
        penc(s)         Primary energy goods               / coa, cru, ngd, per /
* Non-electricity subset
        nele(s)      Non-electricity sectors               / apa, cba, frs, grn, oca,
                                                             osa, pfb, vna, coa, cru,
                                                             cem, con, min, pfm, trm,
                                                             fbm, ngd, bom, bos, cpm,
                                                             tec, wpm, per, fin, hlt,
                                                             pub, rtl, tel, trn /
;

* ---------------------------------------------------------------------------- *
*       Read aggregate mapping set:
* ---------------------------------------------------------------------------- *

$include %load_map%

$include %load_map2%

        DISPLAY s_a, zones, econr, econz_map, zones_map, econr_map, r_psm_zones, r_psm_econr;

* ---------------------------------------------------------------------------- *
*       Install parameters for model specification:
* ---------------------------------------------------------------------------- *

SCALAR
        calib                   Calibration flag - 0 == do not use constraint on consumption
        obj_solution            Benchmark calibration solution
;

PARAMETER
        vdifm(r,g,s)            Total intermediate demand
        evo(r,i,f)              Factor endowment by institution
        vfm_k(r,s)              Capital payments by institution
        vfm_l(r,s)              Wage payments to household by region
        evo_k(r,h)              Capital endowment to household by region
        evo_l(r,h)              Labor endowment by household
        va(r,s)                 Armington supply including imports
        vn(g)                   Intra-national trade
        vinvd(r,g)              Investment demand by commodity
        vinvh(r,h)              Investment demand by household
        incadj(r,h)             Base year net transfer
;

vdifm(r,g,s) = vdfm(r,g,s) + sum(trd, vifm(r,g,trd,s));
evo(r,i,f)   = sum(t, evom(r,f,i,t));
va(r,s)      = sum(trd,vim(r,s,trd)) +  vdmi(r,s);
vn(g)        = sum(r, vxm(r,g,"dtrd"));

        DISPLAY vdmi, vxm, vim, vdifm, evo, va, vn;

* ---------------------------------------------------------------------------- *
*       Impute investment demand to households:
* ---------------------------------------------------------------------------- *

vinvd(r,g) = va(r,g) - (sum(h,vdpm(r,g,h))+sum((trd,h),vipm(r,g,trd,h)))
                     - (sum(pub,vdgm(r,g,pub))+sum((trd,pub),vigm(r,g,trd,pub)))
                     - sum(s, vdifm(r,g,s));

        DISPLAY vinvd;

* ---------------------------------------------------------------------------- *
*       Impute wage payments, capital payments:
* ---------------------------------------------------------------------------- *

vfm(r,"empl",s) = max(vdmi(r,s)+sum(trd,vxm(r,s,trd))
                  - (vfm(r,"prop",s) + vfm(r,"othp",s)) - vfm(r,"btax",s)
                  - sum(g, vdifm(r,g,s)), 0);

vfm(r,"prop",s) = max(vdmi(r,s)+sum(trd,vxm(r,s,trd))
                  - (vfm(r,"empl",s) + vfm(r,"othp",s)) - vfm(r,"btax",s)
                  - sum(g, vdifm(r,g,s)), 0);

vfm(r,"othp",s) = max(vdmi(r,s)+sum(trd,vxm(r,s,trd))
                  - (vfm(r,"empl",s) + vfm(r,"prop",s)) - vfm(r,"btax",s)
                  - sum(g, vdifm(r,g,s)), 0);

vfm(r,"btax",s) = vdmi(r,s)+sum(trd,vxm(r,s,trd))
                  - (vfm(r,"prop",s) + vfm(r,"othp",s)) - vfm(r,"empl",s)
                  - sum(g, vdifm(r,g,s));

        DISPLAY vfm;

* ---------------------------------------------------------------------------- *
*       Scale labor endowments to match imputed wage payments:
* ---------------------------------------------------------------------------- *

evo(r,h,"empl") = evo(r,h,"empl")/sum(hh, evo(r,hh,"empl"))
                  * sum(s, vfm(r,"empl",s));

evo(r,h,"prop") = evo(r,h,"prop")/sum(hh,evo(r,hh,"prop"))
                  * sum(s, vfm(r,"prop",s));

evo(r,h,"othp") = evo(r,h,"othp")/sum(hh,evo(r,hh,"othp"))
                  * sum(s, vfm(r,"othp",s));

        DISPLAY evo;

* ---------------------------------------------------------------------------- *
*        Impute wage payments, capital payments for stock values:
* ---------------------------------------------------------------------------- *

vfm_k(r,s) = vfm(r,"othp",s) + (1)*vfm(r,"prop",s);

vfm_l(r,s) = vfm(r,"empl",s) + (0)*vfm(r,"prop",s);

evo_k(r,h) = evo(r,h,"othp") + (1)*evo(r,h,"prop");

evo_l(r,h) = evo(r,h,"empl") + (0)*evo(r,h,"prop");

        DISPLAY vfm_k, vfm_l, evo_k, evo_l;

* ---------------------------------------------------------------------------- *
*       Impute aggregate investment and investment demand by household:
* ---------------------------------------------------------------------------- *

vinv(r) = sum(g, vinvd(r,g));
vinvh(r,h)$sum(hh, max(0, sum(t,trnsfer(r,"inv",t,hh) - trnsfer(r,hh,t,"inv"))))
                = max(0, sum(t, trnsfer(r,"inv",t,h) - trnsfer(r,h,t,"inv"))) /
                  sum(hh, max(0, sum(t, trnsfer(r,"inv",t,hh)
                  - trnsfer(r,hh,t,"inv"))))* vinv(r);

        DISPLAY vinv, vinvh;

* ---------------------------------------------------------------------------- *
*       Adjust household income:
* ---------------------------------------------------------------------------- *

incadj(r,h) = sum(g,vdpm(r,g,h)+sum(trd,vipm(r,g,trd,h)))
              + vinvh(r,h) - sum(f,evo(r,h,f));

        DISPLAY incadj;

* ---------------------------------------------------------------------------- *
*        Define elasticity parameters:
* ---------------------------------------------------------------------------- *
PARAMETER
* Sectoral elasticities
        es_kl(s)              Elasticity of substitution - value-added
        es_kle(s)             Elasticity of substitution - value-added-energy bundle
        es_klem(s)            Elasticity of substitution - value-added-energy-materials bundle
        es_mat                Elasticity of substitution - materials bundle
        es_ene(s)             Elasticity of substitution - electricity and primary energy bundle
        es_en(s)              Elasticity of substitution - primary energy bundle
        es_rklem(s)           Elasticity of substitution - r-klem bundle in resource extraction sectors
* Household elasticities
        es_c                  Elasticity of substitution - consumption bundle
        es_cen                Elasticity of substitution - consumption electricity and energy bundle
        es_cene               Elasticity of substitution - consumption primary energy bundle
        es_cm                 Elasticity of substitution - consumption materials bundle
* Armington  elasticities
        es_dn(s)              Elasticity of substitution - domestic and imported goods
        es_nf(s)              Elasticity of substitution - domestic and national
* Elasticity of Transformation
        et_dx                 Elasticity of transformation - domestic and exported production
* Resource Elasticities
        res_share(s)          Fixed resource extraction share
        es_res_lr(s)          Elasticity of transformation - long-run resource extraction elasticity
;

$GDXIN %load_es%
$load es_kl es_kle es_klem es_mat es_ene es_en
$load es_c es_cm es_cen es_cene
$load es_dn es_nf et_dx
*$load res_share es_res_lr

        DISPLAY es_kl, es_kle, es_klem, es_mat, es_ene, es_en,
                es_c, es_cm, es_cen, es_cene,
                es_dn, es_nf, et_dx ;

* ---------------------------------------------------------------------------- *
*        Define inital model parameters for simulations:
* ---------------------------------------------------------------------------- *

PARAMETER
        vdmi0(r,s)         Initial domestic output including institutional imake
        vdgm0(r,s,pub)     Initial domestic public demand
        vigm0(r,s,trd,pub) Initial imported public demand
        vipm0(r,s,trd,h)   Initial imported consumption demand
        vdifm0(r,g,s)      Initial total intermediate demand
        vfm0(r,fa,s)       Initial factor demand
        vdmi0(r,s)         Initial domestic output including institutional imake
        evo0(r,i,f)        Initial factor endowment by institution
        vfm_k0(r,s)        Initial capital payments by institution
        vfm_l0(r,s)        Initial wage payments to household by region
        evo_k0(r,h)        Initial capital endowment to household by region
        evo_l0(r,h)        Initial labor endowment by household
        va0(r,s)           Initial armington supply including imports
        vn0(g)             Initial intra-national trade
        vxm0(r,g,trd)      Initial national and international exports
        vinv0(r)           Initial aggregate investment
        vinvd0(r,g)        Initial investment demand by commodity
        vinvh0(r,h)        Initial investment demand by household
        incadj0(r,h)       Initial base year net transfer
;

vdmi0(r,s)         = vdmi(r,s);
vdgm0(r,s,pub)     = vdgm(r,s,pub);
vigm0(r,s,trd,pub) = vigm(r,s,trd,pub);
vipm0(r,s,trd,h)   = vipm(r,s,trd,h);
vdifm0(r,g,s)      = vdifm(r,g,s);
vfm0(r,fa,s)       = vfm(r,fa,s);
evo0(r,i,f)        = evo(r,i,f);
vfm_k0(r,s)        = vfm_k(r,s);
vfm_l0(r,s)        = vfm_l(r,s);
evo_k0(r,h)        = evo_k(r,h);
evo_l0(r,h)        = evo_l(r,h);
va0(r,s)           = va(r,s);
vn0(g)             = vn(g) ;
vinv0(r)           = vinv(r);
vinvd0(r,g)        = vinvd(r,g);
vinvh0(r,h)        = vinvh(r,h);
incadj0(r,h)       = incadj(r,h);
vxm0(r,g,trd)      = vxm(r,g,trd);

* ---------------------------------------------------------------------------- *
*       Formulate SOE US regional model in GAMS/MPSGE based on hese statistics
*       and verify that the resulting data represents an equilibrium:
* ---------------------------------------------------------------------------- *
$ontext

$model:soe

*       Debug options
*$funlog:.true
*$echop:.true
*$funlog:.true. is a switch to generate a det_dxiled listing of function evaluations for all production sectors and consumers.
*$echop:.true. is a switch for returning all or part of the scalar MPSGE file to the solver status file. (enter the GAMS statement "OPTION SYSOUT=ON;" prior to solving the model)
*$PEPS:1.0E-15
*$DATECH:.true
*$EULCHK:.true
*$WALCHK:.true

*               Define prices, commodies and consumers
$sectors:
        y(r,s)$(vdmi(r,s)+sum(trd,vxm(r,s,trd))) ! Sectoral production
        a(r,s)$va(r,s)                           ! Armington aggregation
        c(r,h)                                   ! Consumption by household
        gov(r,pub)                               ! Public output
        inv(r)                                   ! Investment

$commodities:
        p(r,s)$(vdmi(r,s))                       ! Sectoral output prices
        pc(r,h)                                  ! Consumption by household
        pa(r,s)$va(r,s)                          ! Armington aggregate prices
        pn(s)$vn(s)                              ! Intra-national trade price
        pinv(r)                                  ! New investment
        pgov(r,pub)                              ! Public output
        pf_k(r)$(sum(s,vfm_k(r,s)))              ! Factor capital prices
        pf_l(r)$(sum(s,vfm_l(r,s)))              ! Factor labor prices
        pfx                                      ! Foreign exchange
        ptax(r)$(sum(s,vfm(r,"btax",s)))         ! Business taxes

$consumers:
        rh(r,h)                                  ! Representative households
        govt(r,pub)                              ! Government (different levels)
        taxrev(r)                                ! Tax revenue agent

$auxiliary:
        rhscale(r,h)$(not calib) ! Rationing constraint (scales transfers and investment with households activity levels)

*               Define production function for agricultural sectors
$prod:y(r,agrs)$(vdmi(r,agrs)+sum(trd,vxm(r,agrs,trd))) t:et_dx(agrs) s:es_klem(agrs)
+                                            mat:es_mat vae:es_kle(agrs)
+                                            en(vae):es_ene(agrs) va(vae):es_kl(agrs)
+                                            pe(en):es_en(agrs)
        o:p(r,agrs)$(vdmi(r,agrs))   q:vdmi(r,agrs)
        o:pfx                        q:vxm(r,agrs,"ftrd")
        o:pn(agrs)                   q:vxm(r,agrs,"dtrd")
        i:pa(r,g)                    q:((psmtfp(r,g,agrs)*dre_prod(r,g,agrs)*wdrgtfp(r,agrs)*tfp(r,agrs))*vdifm(r,g,agrs))  mat:$matg(g) en:$eleg(g) pe:$peng(g)
        i:pf_k(r)$vfm_k(r,agrs)      q:((psmtfp_fa(r,agrs)*wdrgtfp(r,agrs)*tfp(r,agrs))*vfm_k(r,agrs))                      va:
        i:pf_l(r)$vfm_l(r,agrs)      q:((psmtfp_fa(r,agrs)*wdrgtfp(r,agrs)*tfp(r,agrs))*vfm_l(r,agrs))                      va:
        i:ptax(r)                    q:vfm(r,"btax",agrs)

*               Define production function for material sectors (construction, electricity, manufacturing, petrolium refineries services, and transportation)
$prod:y(r,mats)$(vdmi(r,mats)+sum(trd,vxm(r,mats,trd))) t:et_dx(mats) s:es_klem(mats)
+                                            mat:es_mat vae:es_kle(mats)
+                                            en(vae):es_ene(mats) va(vae):es_kl(mats)
+                                            pe(en):es_en(mats)
        o:p(r,mats)$(vdmi(r,mats))   q:vdmi(r,mats)
        o:pfx                        q:vxm(r,mats,"ftrd")
        o:pn(mats)                   q:vxm(r,mats,"dtrd")
        i:pa(r,g)                    q:((psmtfp(r,g,mats)*dre_prod(r,g,mats)*tfp(r,mats))*vdifm(r,g,mats))  mat:$matg(g) en:$eleg(g) pe:$peng(g)
        i:pf_k(r)$vfm_k(r,mats)      q:((psmtfp_fa(r,mats)*tfp(r,mats))*vfm_k(r,mats))                      va:
        i:pf_l(r)$vfm_l(r,mats)      q:((psmtfp_fa(r,mats)*tfp(r,mats))*vfm_l(r,mats))                      va:
        i:ptax(r)                    q:vfm(r,"btax",mats)

*               Define production function for resource extraction sectors sectors
$prod:y(r,rexs)$(vdmi(r,rexs)+sum(trd,vxm(r,rexs,trd))) t:et_dx(rexs) s:es_klem(rexs)
+                                            mat:es_mat  vae:es_kle(rexs)
+                                            en(vae):es_ene(rexs)  va(vae):es_kl(rexs)
+                                            pe(en):es_en(rexs)
        o:p(r,rexs)$(vdmi(r,rexs))   q:vdmi(r,rexs)
        o:pfx                        q:vxm(r,rexs,"ftrd")
        o:pn(rexs)                   q:vxm(r,rexs,"dtrd")
        i:pa(r,g)                    q:((psmtfp(r,g,rexs)*dre_prod(r,g,rexs)*tfp(r,rexs))*vdifm(r,g,rexs))  mat:$matg(g) en:$eleg(g) pe:$peng(g)
        i:pf_k(r)$vfm_k(r,rexs)      q:((psmtfp_fa(r,rexs)*tfp(r,rexs))*vfm_k(r,rexs))                      va:
        i:pf_l(r)$vfm_l(r,rexs)      q:((psmtfp_fa(r,rexs)*tfp(r,rexs))*vfm_l(r,rexs))                      va:
        i:ptax(r)                    q:vfm(r,"btax",rexs)

*               Define Armington trade production functions
$prod:a(r,s)$va(r,s)  s:es_nf(s)  m:es_dn(s)
        o:pa(r,s)                    q:va(r,s)
        i:p(r,s)$(vdmi(r,s))         q:vdmi(r,s)                m:
        i:pfx                        q:vim(r,s,"ftrd")
        i:pn(s)                      q:vim(r,s,"dtrd")          m:

*               Define investment production function
$prod:inv(r)
        o:pinv(r)                    q:vinv(r)
        i:pa(r,s)                    q:vinvd(r,s)

*               Define government production function
$prod:gov(r,pub)
        o:pgov(r,pub)                q:vgm(r,pub)
        i:pa(r,s)                    q:(vdgm(r,s,pub)+sum(trd,vigm(r,s,trd,pub)))

*               Define consumption demand production function
*$prod:c(r,h) s:1
*        o:pc(r,h)                    q:vpm(r,h)
*        i:pa(r,s)                    q:(dre_hh(r,s,h)*vdpm(r,s,h)+sum(trd,vipm(r,s,trd,h)))
*               Define consumption demand production function
$prod:c(r,h) s:es_c
+                mat:es_cm ene:es_cene pe(ene):es_cen
        o:pc(r,h)                    q:vpm(r,h)
        i:pa(r,s)                    q:(dre_hh(r,s,h)*vdpm(r,s,h)+sum(trd,vipm(r,s,trd,h))) mat:$matc(s) ene:$elec(s) pe:$penc(s)

*               Define household utility function
$demand:rh(r,h)
        d:pc(r,h)                    q:vpm(r,h)
        e:pf_k(r)                    q:evo_k(r,h)
        e:pf_l(r)                    q:evo_l(r,h)
        e:pfx                        q:incadj(r,h)              r:rhscale(r,h)$(not calib)
        e:pinv(r)                    q:(-vinvh(r,h))            r:rhscale(r,h)$(not calib)

*               Define government demand function
$demand:govt(r,pub)
        d:pgov(r,pub)                q:vgm(r,pub)
        e:pfx                        q:vgm(r,pub)

*               Define tax revenue function
$demand:taxrev(r)
        d:pfx
        e:ptax(r)                    q:(sum(s,vfm(r,"btax",s)))

*               Defind household constraint on consumption through income and investment
$constraint:rhscale(r,h)$(not calib)

        rhscale(r,h)    =e=          c(r,h);

*               Define model output parameters
$report:
*  Production domestic, export intranational, export international
        v:yd(r,s)           o:p(r,s)        prod:y(r,s)
        v:yxf(r,s)          o:pfx           prod:y(r,s)
        v:yxd(r,s)          o:pn(s)         prod:y(r,s)
        v:id(r,g,s)         i:pa(r,g)       prod:y(r,s)
        v:kd(r,s)           i:pf_k(r)       prod:y(r,s)
        v:ld(r,s)           i:pf_l(r)       prod:y(r,s)
        v:btd(r,f,s)        i:ptax(r)       prod:y(r,s)
* GRP consumption, government. investment demand
        v:cd(r,h)           o:pc(r,h)       prod:c(r,h)
        v:govd(r,pub)       d:pgov(r,pub)   demand:govt(r,pub)
        v:invd(r)           o:pinv(r)       prod:inv(r)
* Armington, imports (foreign, imports (domestic)
        v:aa(r,s)           o:pa(r,s)       prod:a(r,s)
        v:mf(r,s)           i:pfx           prod:a(r,s)
        v:md(r,s)           i:pn(s)         prod:a(r,s)
*
        v:evolp(r,h)        d:pf_l(r)       demand:rh(r,h)
        v:evokp(r,h)        d:pf_k(r)       demand:rh(r,h)
        v:hinv(r,h)         d:pinv(r)       demand:rh(r,h)
        v:ca(r,s,h)         i:pa(r,s)       prod:c(r,h)
$offtext

* ---------------------------------------------------------------------------- *
*       Run model calibtration:
* ---------------------------------------------------------------------------- *

*       Call MPSGE program:

$sysinclude mpsgeset soe

*       Choose a numeraire:

PFX.fx = 1;

*       Define MPSGE solve parameters

soe.workspace = 40;
soe.iterlim = 0;
calib = 1;
$include SOE.GEN
solve soe using mcp;

*       Exit if model does not calibrate

obj_solution  = abs(soe.objval);

        DISPLAY obj_solution;

abort$(abs(soe.objval) gt 8e-5) "***Model does not calibrate***";

* ---------------------------------------------------------------------------- *
*       Define sets on  policy/shock simulations:
* ---------------------------------------------------------------------------- *

SET
* Main recursive sets
       tp               Year 2010-2100             / 2010*2100 /
       sim              simulation variables       / bau             "Business-As-Usual (no simulation shocks)",
                                                     cntf_wdrg       "WDRG Impact Simulation",
                                                     cntf_psm_nosul  "PSM Electricity Impacts Simulation (without other productivity impacts)",
                                                     cntf_psm_full   "PSM Electricity and Productivity Impacts Simulation",
                                                     cntf_sts        "Demand Response Emulation (STS) Impacts Simulation",
                                                     psm_nosul_sts   "PSM Electricity and Demand Response (STS) Impacts Simulation",
                                                     psm_full_sts    "PSM Electricity/Productivity and Demand Response (STS) Impacts Simulation",
                                                     cntf_td         "Demand Response Emulation (TD) Impacts Simulation",
                                                     psm_nosul_td    "PSM Electricity and Demand Response (TD) Impacts Simulation",
                                                     psm_full_td     "PSM Electricity/Productivity and Demand Response (TD) Impacts Simulation" /
;

ALIAS (sim, ssim) ;

* ---------------------------------------------------------------------------- *
*       Define sets for dynamic recursive loop with policy/shock simulations:
* ---------------------------------------------------------------------------- *

SET
* Time period subset
       tp_2010(tp)              Time period for the running the static model
                                / 2010 /
       tp_2015(tp)              Time periods 2010-2015 in five year steps (9 tp)
                                / 2010, 2015 /
       tp_2020(tp)              Time periods 2010-2020 in five year steps (9 tp)
                                / 2010, 2015, 2020 /
       tp_2025(tp)              Time periods 2010-2025 in five year steps (9 tp)
                                / 2010, 2015, 2020, 2025 /
       tp_2030(tp)              Time periods 2010-2030 in five year steps (9 tp)
                                / 2010, 2015, 2020, 2025, 2030 /
       tp_2035(tp)              Time periods 2010-2035 in five year steps (9 tp)
                                / 2010, 2015, 2020, 2025, 2030,
                                  2035 /
       tp_2040(tp)              Time periods 2010-2040 in five year steps (9 tp)
                                / 2010, 2015, 2020, 2025, 2030,
                                  2035, 2040 /
       tp_2045(tp)              Time periods 2010-2045 in five year steps (9 tp)
                                / 2010, 2015, 2020, 2025, 2030,
                                  2035, 2040, 2045 /
       tp_2050(tp)              Time periods 2010-2050 in five year steps (9 tp)
                                / 2010, 2015, 2020, 2025, 2030,
                                  2035, 2040, 2045, 2050 /

       tp_2100(tp)              Time periods 2010-2100 in five year steps (19 tp)
                                / 2010, 2015, 2020, 2025, 2030,
                                  2035, 2040, 2045, 2050, 2055,
                                  2060, 2065, 2070, 2075, 2080,
                                  2085, 2090, 2095, 2100 /

* Coupled Model subsets
       sim_bau(sim)             Simulation variables for business-as-usual
                                / bau /
       sim_wdrg(sim)             Simulation variables for WDRG coupling
                                / cntf_wdrg /
       sim_psm(sim)             Simulation variables for PSM (non-productivity) coupling
                                / cntf_psm_nosul /
       sim_psmfull(sim)         Simulation variables for PSM (w\productivity) coupling
                                / cntf_psm_full /
       sim_drm_sts(sim)         Simulation variables for DRM coupling
                                / cntf_sts /
       sim_psmdrm(sim)          Simulation variables for Full PSM and DRM coupling
                                / psm_nosul_sts /
       sim_psmfulldrm(sim)      Simulation variables for Full PSM and DRM coupling
                                / psm_full_sts /
* Dynamic subsets
       simtp(tp)                Time periods of simulation
       policy(sim)              Policy simulations of multiple coupled models
       policy_bau(sim)          BAU with single policy shock
       forloop(sim)             Policy simulations of recursive loop
;

* ---------------------------------------------------------------------------- *
*        Define time frame years based on switch selection:
* ---------------------------------------------------------------------------- *

if (switch_time = 2010,
        simtp(tp)$tp_2010(tp) = Yes;
elseif (switch_time = 2015),
        simtp(tp)$tp_2015(tp) = Yes;
elseif (switch_time = 2015),
        simtp(tp)$tp_2020(tp) = Yes;
elseif (switch_time = 2020),
        simtp(tp)$tp_2025(tp) = Yes;
elseif (switch_time = 2025),
        simtp(tp)$tp_2030(tp) = Yes;
elseif (switch_time = 2030),
        simtp(tp)$tp_2035(tp) = Yes;
elseif (switch_time = 2035),
        simtp(tp)$tp_2040(tp) = Yes;
elseif (switch_time = 2040),
        simtp(tp)$tp_2040(tp) = Yes;
elseif (switch_time = 2045),
        simtp(tp)$tp_2045(tp) = Yes;
elseif (switch_time = 2050),
        simtp(tp)$tp_2050(tp) = Yes;
elseif (switch_time = 2100),
        simtp(tp)$tp_2100(tp) = Yes;
);
        DISPLAY simtp;

* ---------------------------------------------------------------------------- *
*        Define time steps between time periods and over the whole period:
* ---------------------------------------------------------------------------- *

SET
        tfirst(tp)    First simualtion year
        tlast(tp)     Last simulation year
;

PARAMETER
        tsval(tp)     Simulation time periods as values
        ts(tp)        Number of year between each time period
;

SCALAR
        tts           Total number of years over the whole simulation time frame
;

* Relax ordered set for dynamic sets (keep off to use simtp)
$offOrder

* Define first and last time period date
tfirst(simtp) = yes$(ord(simtp) = 1);
tlast(simtp) = yes$(ord(simtp) = card(simtp));

* Define time step
tsval(simtp) = simtp.val;
*ts(simtp)$(not simtp.last) = tsval(simtp+1)-tsval(simtp);
ts(simtp+1) = tsval(simtp+1)-tsval(simtp);
tts = sum(simtp, ts(simtp));

        DISPLAY tfirst, tlast, tsval, ts, tts;

* ---------------------------------------------------------------------------- *
*        Define time policy simulation based on coupled model selection:
* ---------------------------------------------------------------------------- *

if (   (switch_wdrg = 0 and
        switch_psm = 0 and
        switch_psmf = 0 and
        switch_dre = 0),

            policy_bau(sim)$sim_bau(sim)      = YES;
            policy(sim)$sim_bau(sim)          = YES;

elseif (switch_wdrg = 1),

            policy_bau(sim)$sim_bau(sim)      = YES;
            policy(sim)$sim_wdrg(sim)         = YES;

elseif (switch_psm = 1 and
        switch_psmf = 0 and
        switch_dre = 0),

            policy_bau(sim)$sim_bau(sim)      = YES;
            policy(sim)$sim_psm(sim)          = YES;

elseif (switch_psm = 1 and
        switch_psmf = 1 and
        switch_dre = 0),

            policy_bau(sim)$sim_bau(sim)      = YES;
            policy(sim)$sim_psmfull(sim)      = YES;

elseif (switch_psm = 0 and
        switch_psmf = 0 and
        switch_dre = 1),

            policy_bau(sim)$sim_bau(sim)      = YES;
            policy(sim)$sim_drm_sts(sim)      = YES;

elseif (switch_psm = 1 and
        switch_psmf = 0 and
        switch_dre = 1),

            policy_bau(sim)$sim_bau(sim)      = YES;
            policy_bau(sim)$sim_drm_sts(sim)  = YES;
            policy(sim)$sim_psmdrm(sim)       = YES;

elseif (switch_psm = 1 and
        switch_psmf = 1 and
        switch_dre = 1),

            policy_bau(sim)$sim_bau(sim)      = YES;
            policy_bau(sim)$sim_drm_sts(sim)  = YES;
            policy(sim)$sim_psmfulldrm(sim)   = YES;
);

        DISPLAY policy;

* ---------------------------------------------------------------------------- *
*        Define policy simulation for loop:
* ---------------------------------------------------------------------------- *

forloop(sim)$(policy_bau(sim)) = YES;
forloop(sim)$(policy(sim))     = YES;

        DISPLAY forloop;

* ---------------------------------------------------------------------------- *
*        Define parameters for dynamic model:
* ---------------------------------------------------------------------------- *

PARAMETER
* Input Parameter
         us_pop_proj(r,tp)                      U.S population projections based on SSP
         us_pop_proj_cagr(r,tp)                 U.S population projections Growth Rates - region r

* Parameters for recursive loop
         lstock0(r)                             Base year Labor endowment - region r
         kstock0(r)                             Base year Capital endowment - region r
         lstock(r,tp)                           Labor stock - region r
         kstock(r,tp)                           Capital stock - region r
         evoi_k0(r,h)                           Base year capital endowment - region r for household h
         evoi_l0(r,h)                           Base year labor endowment - region r for household h
         kstock_sim(r,tp,sim)                   Capital stock - region r at time period tp under physical shock sim
         lstock_sim(r,tp,sim)                   Labor stock - region r at time period tp under physical shock sim
         drate                                  Depreciation rate
         irate                                  Interest rate
         kprate                                 National capital productivity
         lprate                                 National labor productivity (based on SAGE AEO 2018 -- 0.016)
         lpidx(r,tp)                            Labor productivity index - region r
         poprate(r,tp)                          National population growt rate (compound annualized growth rate)

* Parameters checks from recursive loop
         lpidx_sim(r,tp)                        Labor productivity index by time period
         psmtfp_sim(r,g,s,tp,sim)               PSM's chg in average marginal cost from shock in region r by time period
         psmtfp_fa_sim(r,s,tp,sim)              PSM's factor chg in average marginal cost from shock in region r by time period
         shk_dr_pchg_prod_sim(r,g,s,tp,sim)     Demand response percentage change in per capita load by time period
         shk_dr_pchg_consum_sim(r,s,h,tp,sim)   Demand response percentage change in per capita load by time period
         evok_sim(r,h,tp,sim)                   Capital endowment -region r for household h at time period tp under physical shock sim
         evol_sim(r,h,tp,sim)                   Labor endowment - region r for household h at time period tp under physical shock sim
         incadj_sim(r,h,tp,sim)
         rhscale_sim(r,h,tp,sim)
         evolp_sim(r,h,tp,sim)
         evokp_sim(r,h,tp,sim)
         ca_sim(r,s,h,tp,sim)
         hinv_sim(r,h,tp,sim)
         incadj_cal_sim(r,h,tp,sim)
;

* ---------------------------------------------------------------------------- *
*       Define parameters for dynamic output - level and percentage values:
* ---------------------------------------------------------------------------- *

PARAMETER
* Prices
        py_sim(r,s,tp,sim)           price - sectoral output price - region r for sector s at time period tp under physical shock sim
        pa_sim(r,s,tp,sim)           price - armington aggregate price - region r for sector s at time period tp under physical shock sim
        pc_sim(r,h,tp,sim)           price - consumption price - region r for sector s at time period 2050 under physical shock sim
        pn_sim(s,tp,sim)             price - intranational trade price - region r for sector s at time period 2050 under physical shock sim
        pinv_sim(r,tp,sim)           price - new investment price - region r for sector s at time period 2050 under physical shock sim
        pgov_sim(r,pub,tp,sim)       price - public output price - region r for sector s at time period 2050 under physical shock sim
        pk_sim(r,tp,sim)             price - factor capital price - region r at time period tp under physical shock sim
        pl_sim(r,tp,sim)             price - factor labor price - region r at time period tp under physical shock sim
        pbt_sim(r,tp,sim)            price - business taxes - region r at time period tp under physical shock sim

        py_pchg(r,s,sim,ssim)        percentage change from the baseline - sectoral output price - region r for sector s at time period 2050 under physical shock sim
        pa_pchg(r,s,sim,ssim)        percentage change from the baseline - armington aggregate price - region r for sector s at time period 2050 under physical shock sim
        pc_pchg(r,h,sim,ssim)        percentage change from the baseline - consumption price - region r for household h at time period 2050 under physical shock sim
        pn_pchg(s,sim,ssim)          percentage change from the baseline - intranational trade price - sector s at time period 2050 under physical shock sim
        pinv_pchg(r,sim,ssim)        percentage change from the baseline - new investment price - region r at time period 2050 under physical shock sim
        pgov_pchg(r,pub,sim,ssim)    percentage change from the baseline - public output price - region r at time period 2050 under physical shock sim
        pk_pchg(r,sim,ssim)          percentage change from the baseline - factor capital price - region r at time period 2050 under physical shock sim
        pl_pchg(r,sim,ssim)          percentage change from the baseline - factor labor price - region r at time period 2050 under physical shock sim
        pbt_pchg(r,sim,ssim)         percentage change from the baseline - business taxes - region r at time period 2050 under physical shock sim

* Activity
        yd_sim(r,s,tp,sim)           activity - sectoral production for the domestic market - region r for sector s at time period tp under physical shock sim
        yd_sim_r(r,tp,sim)           activity - sectoral production for the domestic market - region r at time period tp under physical shock sim
        yd_sim_us(tp,sim)            activity - sectoral production for the domestic market - time period tp under physical shock sim
        vyd_sim(r,s,tp,sim)          activity - sectoral production (value) for the domestic market - region r for sector s at time period tp under physical shock sim
        yxd_sim(r,s,tp,sim)          activity - sectoral production for the intranational market - region r for sector s at time period tp under physical shock sim
        yxd_sim_r(r,tp,sim)          activity - sectoral production for the intranational market - region r at time period tp under physical shock sim
        yxd_sim_us(tp,sim)           activity - sectoral production for the intranational market - time period tp under physical shock sim
        yxf_sim(r,s,tp,sim)          activity - sectoral production for the international market - region r for sector s at time period tp under physical shock sim
        yxf_sim_r(r,tp,sim)          activity - sectoral production for the international  market - region r at time period tp under physical shock sim
        yxf_sim_us(tp,sim)           activity - sectoral production for the international market - time period tp under physical shock sim
        yt_sim(r,s,tp,sim)           activity - total sectoral production for all markets - region r for sector s at time period tp under physical shock sim
        vyt_sim(r,s,tp,sim)          activity - total sectoral production (value) all markets - region r for sector s at time period tp under physical shock sim
        yt_sim_r(r,tp,sim)           activity - total sectoral production for all markets - region r at time period tp under physical shock sim
        yt_sim_us(tp,sim)            activity - total sectoral production for all markets - time period tp under physical shock sim
        a_sim(r,s,tp,sim)            activity - armington aggregation - region r for sector s at time period tp under physical shock sim

        yt_share_sa(r,s,tp,sim)      activity - sub-sector share of total production of aggregate sector - region r sectors for sector s at time period tp under physical shock sim
        yt_share_r(r,s,tp,sim)       activity - sub-sector share of total production of national total production - region r sectors for sector s at time period tp under physical shock sim

        yd_pchg(r,s,sim,ssim)        percentage change from the baseline - sectoral production for the domestic market - region r for sector s at time period 2050 under physical shock sim
        yd_pchg_r(r,sim,ssim)        percentage change from the baseline - sectoral production for the domestic market - region r at time period 2050 under physical shock sim
        yd_pchg_us(sim,ssim)         percentage change from the baseline - sectoral production for the domestic market - region r at time period 2050 under physical shock sim
        vyd_pchg(r,s,sim,ssim)       percentage change from the baseline - sectoral production (value) for the domestic market - region r for sector s at time period 2050 under physical shock sim
        yxd_pchg(r,s,sim,ssim)       percentage change from the baseline - sectoral production for the intranational market - region r for sector s at time period 2050 under physical shock sim
        yxd_pchg_r(r,sim,ssim)       percentage change from the baseline - sectoral production for the intranational market - region r at time period 2050 under physical shock sim
        yxd_pchg_us(sim,ssim)        percentage change from the baseline - sectoral production for the intranational market - time period 2050 under physical shock sim
        yxf_pchg(r,s,sim,ssim)       percentage change from the baseline - sectoral production for the international market - region r for sector s at time period 2050 under physical shock sim
        yxf_pchg_r(r,sim,ssim)       percentage change from the baseline - sectoral production for the international market - region r at time period 2050 under physical shock sim
        yxf_pchg_us(sim,ssim)        percentage change from the baseline - sectoral production for the international  market - time period 2050 under physical shock sim
        yt_pchg(r,s,sim,ssim)        percentage change from the baseline - total sectoral production for the domestic market - region r for sector s at time period 2050 under physical shock sim
        yt_pchg_r(r,sim,ssim)        percentage change from the baseline - total sectoral production for the domestic market - region r at time period 2050 under physical shock sim
        yt_pchg_us(sim,ssim)         percentage change from the baseline - total sectoral production for the domestic market - time period 2050 under physical shock sim
        a_pchg(r,s,sim,ssim)         percentage change from the baseline - armington aggregation - region r for sector s at time period 2050 under physical shock sim

* Activity by subset
        yd_sim_sa(r,s_a,tp,sim)      activity - sectoral production for the domestic market - region r for aggregate sector s_a at time period tp under physical shock sim
        yxd_sim_sa(r,s_a,tp,sim)     activity - sectoral production for the intranational market - region r for aggregate sector s_a at time period tp under physical shock sim
        yxf_sim_sa(r,s_a,tp,sim)     activity - sectoral production for the international market - region r for aggregate sector s_a at time period tp under physical shock sim
        yt_sim_sa(r,s_a,tp,sim)      activity - total sectoral production for the domestic market - region r for aggregate sector s_a at time period tp under physical shock sim

        yd_pchg_sa(r,s_a,sim,ssim)   percentage change from the baseline - sectoral production for the domestic market - region r for aggregate sector s_a at time period tp under physical shock sim
        yxd_pchg_sa(r,s_a,sim,ssim)  percentage change from the baseline - sectoral production for the intranational market - region r for aggregate sector s_a at time period tp under physical shock sim
        yxf_pchg_sa(r,s_a,sim,ssim)  percentage change from the baseline - sectoral production for the international market - region r for aggregate sector s_a at time period tp under physical shock sim
        yt_pchg_sa(r,s_a,sim,ssim)   percentage change from the baseline - total sectoral production for the domestic market - region r for aggregate sector s_a at time period tp under physical shock sim

* Demand - factor
        kd_sim(r,s,tp,sim)           demand - capital demand - region r for sector s at time period tp under physical shock sim
        kd_sim_r(r,tp,sim)           demand - payments to capital - region r at time period tp under physical shock sim
        kd_sim_us(tp,sim)            demand - payments to capital - time period tp under physical shock sim
        ld_sim(r,s,tp,sim)           demand - payments to labor (wages) - region r for sector s at time period tp under physical shock sim
        ld_sim_r(r,tp,sim)           demand - payments to labor (wages) - region r at time period tp under physical shock sim
        ld_sim_us(tp,sim)            demand - payments to labor (wages) - time period tp under physical shock sim
        btd_sim(r,s,tp,sim)          demand - payments to (business) taxes - region r for sector s at time period tp under physical shock sim
        btd_sim_r(r,tp,sim)          demand - payments to (business) taxes - region r at time period tp under physical shock sim
        btd_sim_us(tp,sim)           demand - payments to (business) taxes - time period tp under physical shock sim

        kd_pchg(r,s,sim,ssim)        percentage change from the baseline - payments to capital - region r for sector s at time period 2050 under physical shock sim
        kd_pchg_r(r,sim,ssim)        percentage change from the baseline - payments to capital - region r at time period 2050 under physical shock sim
        kd_pchg_us(sim,ssim)         percentage change from the baseline - payments to capital - time period 2050 under physical shock sim
        ld_pchg(r,s,sim,ssim)        percentage change from the baseline - payments to labor (wages) - region r for sector s at time period 2050 under physical shock sim
        ld_pchg_r(r,sim,ssim)        percentage change from the baseline - payments to labor (wages) - region r at time period 2050 under physical shock sim
        ld_pchg_us(sim,ssim)         percentage change from the baseline - payments to labor (wages) - time period 2050 under physical shock sim
        btd_pchg(r,s,sim,ssim)       percentage change from the baseline - payments to (business) taxes - region r for sector s at time period 2050 under physical shock sim
        btd_pchg_r(r,sim,ssim)       percentage change from the baseline - payments to (business) taxes - region r at time period 2050 under physical shock sim
        btd_pchg_us(sim,ssim)        percentage change from the baseline - payments to (business) taxes - time period 2050 under physical shock sim
* Demand - intermediate
        id_sim(r,g,s,tp,sim)         demand - intermediate demand - region r for good g going to sector s at time period tp under physical shock sim
        id_sim_g(g,tp,sim)           demand - intermediate demand -good g at time period tp under physical shock sim
        id_sim_r(r,tp,sim)           demand - intermediate demand - region r at time period tp under physical shock sim
        id_sim_us(tp,sim)            demand - intermediate demand - time period tp under physical shock sim

        id_pchg(r,g,s,sim,ssim)      percentage change from the baseline - intermediate demand - region r for good g going to sector s at time period 2050 under physical shock sim
        id_pchg_g(g,sim,ssim)        percentage change from the baseline - intermediate demand - good g at time period 2050 under physical shock sim
        id_pchg_r(r,sim,ssim)        percentage change from the baseline - intermediate demand - region r at time period tp under physical shock sim
        id_pchg_us(sim,ssim)         percentage change from the baseline - intermediate demand - time period 2050 under physical shock sim
* Demand - export and import
        md_sim(r,s,tp,sim)           demand - imports from intranational trade - region r for sector s at time period tp under physical shock sim
        md_sim_r(r,tp,sim)           demand - imports from intranational trade - region r at time period tp under physical shock sim
        md_sim_us(tp,sim)            demand - imports from intranational trade - time period tp under physical shock sim
        mf_sim(r,s,tp,sim)           demand - imports from international trade - region r for sector s at time period tp under physical shock sim
        mf_sim_r(r,tp,sim)           demand - imports from international trade - region r at time period tp under physical shock sim
        mf_sim_us(tp,sim)            demand - imports from international trade - time period tp under physical shock sim

        md_pchg(r,s,sim,ssim)        percentage change from the baseline - imports from intranational trade - region r for sector s at time period 2050 under physical shock sim
        md_pchg_r(r,sim,ssim)        percentage change from the baseline - imports from intranational trade - region r at time period 2050 under physical shock sim
        md_pchg_us(sim,ssim)         percentage change from the baseline - imports from intranational trade - time period 2050 under physical shock sim
        mf_pchg(r,s,sim,ssim)        percentage change from the baseline - imports from international trade - region r for sector s at time period 2050 under physical shock sim
        mf_pchg_r(r,sim,ssim)        percentage change from the baseline - imports frominternational trade - region r at time period 2050 under physical shock sim
        mf_pchg_us(sim,ssim)         percentage change from the baseline - imports from international trade - time period 2050 under physical shock sim
* Demand - consumption, investiment, and government
        cd_sim(r,h,tp,sim)           demand - Consumption - region r for household h at time period tp under physical shock sim
        cd_sim_h(h,tp,sim)           demand - Consumption - household h at time period tp under physical shock sim
        cd_sim_r(r,tp,sim)           demand - Consumption - region r at time period tp under physical shock sim
        cd_sim_us(tp,sim)            demand - Consumption - time period tp under physical shock sim
        i_sim_r(r,tp,sim)            demand - investment - region r at time period tp under physical shock sim
        i_sim_us(tp,sim)             demand - investment - time period tp under physical shock sim
        g_sim(r,pub,tp,sim)          demand - public output - region r for sector s at time period tp under physical shock sim
        g_sim_r(r,tp,sim)            demand - public output - region r at time period tp under physical shock sim
        g_sim_us(tp,sim)             demand - public output - time period tp under physical shock sim

        cd_pchg(r,h,sim,ssim)        percentage change from the baseline - consumption - region r for household h at time period tp under physical shock sim
        cd_pchg_h(h,sim,ssim)        percentage change from the baseline - consumption - household h at time period tp under physical shock sim
        cd_pchg_r(r,sim,ssim)        percentage change from the baseline - consumption - region r at time period tp under physical shock sim
        cd_pchg_us(sim,ssim)         percentage change from the baseline - consumption - time period tp under physical shock sim
        i_pchg_r(r,sim,ssim)         percentage change from the baseline - investment - region r at time period 2050 under physical shock sim
        i_pchg_us(sim,ssim)          percentage change from the baseline - investment - time period 2050 under physical shock sim
        g_pchg(r,pub,sim,ssim)       percentage change from the baseline - public output - region r for sector s at time period 2050 under physical shock sim
        g_pchg_r(r,sim,ssim)         percentage change from the baseline - public output - region r at time period 2050 under physical shock sim
        g_pchg_us(sim,ssim)          percentage change from the baseline - public output - time period 2050 under physical shock sim
* Tax revenue
        trev_sim_r(r,tp,sim)         additional - tax revenue - region r at time period tp under physical shock sim
        trev_sim_us(tp,sim)          additional - tax revenue - time period tp under physical shock sim

        trev_pchg_r(r,sim,ssim)      percentage change from the baseline - tax revenue - region r at time period 2050 under physical shock sim
        trev_pchg_us(sim,ssim)       percentage change from the baseline - tax revenue - time period 2050 under physical shock sim
* Gross regional product
        grpfc_sim_r(r,tp,sim)        additional - gross Regional Product (final demand) - region r at time period tp under physical shock sim
        grpfc_sim_us(tp,sim)         additional - gross domestic product (final demand) - time period tp under physical shock sim
        grpva_sim_r(r,tp,sim)        additional - gross regional Product (value added) - region r at time period tp under physical shock sim
        grpva_sim_us(tp,sim)         additional - gross domestic poduct (value added) - time period tp under physical shock sim

        grpfc_pchg_r(r,sim,ssim)     percentage change from the baseline - gross regional Product (final demand) - region r at time period 2050 under physical shock sim
        grpfc_pchg_us(sim,ssim)      percentage change from the baseline - gross domestic product (final demand) - physical shock sim
        grpva_pchg_r(r,sim,ssim)     percentage change from the baseline - gross regional Product (value added) - region r at time period 2050 under physical shock sim
        grpva_pchg_us(sim,ssim)      percentage change from the baseline - gross domestic poduct (value added) - physical shock sim
;

* ---------------------------------------------------------------------------- *
*       Define parameters for dynamic output - Compound annualized growth rates:
* ---------------------------------------------------------------------------- *

PARAMETER
* CAGR
        yd_cagr(r,s,tp,sim)          between period cagr - sectoral production for the domestic market - region r for sector s at time period tp under physical shock sim
        yd_cagr_r(r,tp,sim)          between period cagr -  sectoral production for the domestic market - region r at time period tp under physical shock sim
        yd_cagr_us(tp,sim)           between period cagr - sectoral production for the domestic market - time period tp under physical shock sim
        yxd_cagr(r,s,tp,sim)         between period cagr - sectoral production for the intranational market - region r for sector s at time period tp under physical shock sim
        yxd_cagr_r(r,tp,sim)         between period cagr - sectoral production for the intranational market - region r at time period tp under physical shock sim
        yxd_cagr_us(tp,sim)          between period cagr - sectoral production for the intranational market - time period tp under physical shock sim
        yxf_cagr(r,s,tp,sim)         between period cagr - sectoral production for the international market - region r for sector s at time period tp under physical shock sim
        yxf_cagr_r(r,tp,sim)         between period cagr - sectoral production for the international  market - region r at time period tp under physical shock sim
        yxf_cagr_us(tp,sim)          between period cagr - sectoral production for the international market - time period tp under physical shock sim
        yt_cagr(r,s,tp,sim)          between period cagr - total sectoral production for the domestic market - region r for sector s at time period tp under physical shock sim
        yt_cagr_r(r,tp,sim)          between period cagr - total sectoral production for the domestic market - region r at time period tp under physical shock sim
        yt_cagr_us(tp,sim)           between period cagr - total sectoral production for the domestic market - time period tp under physical shock sim

        yd_cagr_sa(r,s_a,tp,sim)     between period cagr - sectoral production for the domestic market - region r for aggregate sector s_a at time period tp under physical shock sim
        yxd_cagr_sa(r,s_a,tp,sim)    between period cagr - sectoral production for the intranational market - region r for aggregate sector s_a at time period tp under physical shock sim
        yxf_cagr_sa(r,s_a,tp,sim)    between period cagr - sectoral production for the international market - region r for aggregate sector s_a at time period tp under physical shock sim
        yt_cagr_sa(r,s_a,tp,sim)     between period cagr - total sectoral production for the domestic market - region r for aggregate sector s_a at time period tp under physical shock sim

        grpfc_cagr_r(r,tp,sim)       between period cagr -  gross regional product (final demand) - region r at time period tp under physical shock sim
        grpfc_cagr_us(tp,sim)        between period cagr - gross domestic product (final demand) - time period tp under physical shock sim
        grpva_cagr_r(r,tp,sim)       between period cagr - gross regional product (value added) - region r at time period tp under physical shock sim
        grpva_cagr_us(tp,sim)        between period cagr - gross domestic poduct (value added) - time period tp under physical shock sim

* Total CAGR
        yd_tcagr(r,s,sim)            cagr over the enitre time frame - sectoral production for the domestic market - region r for sector s under physical shock sim
        yd_tcagr_r(r,sim)            cagr over the enitre time frame - sectoral production for the domestic market - region r under physical shock sim
        yd_tcagr_us(sim)             cagr over the enitre time frame - sectoral production for the domestic market - time period tp under physical shock sim
        yxd_tcagr(r,s,sim)           cagr over the enitre time frame - sectoral production for the intranational market - region r for sector s under physical shock sim
        yxd_tcagr_r(r,sim)           cagr over the enitre time frame - sectoral production for the intranational market - region r under physical shock sim
        yxd_tcagr_us(sim)            cagr over the enitre time frame - sectoral production for the intranational market under physical shock sim
        yxf_tcagr(r,s,sim)           cagr over the enitre time frame - sectoral production for the international market - region r for sector s under physical shock sim
        yxf_tcagr_r(r,sim)           cagr over the enitre time frame - sectoral production for the international  market - region r under physical shock sim
        yxf_tcagr_us(sim)            cagr over the enitre time frame - sectoral production for the international market - time period tp under physical shock sim
        yt_tcagr(r,s,sim)            cagr over the enitre time frame - total sectoral production for the domestic market - region r for sector s under physical shock sim
        yt_tcagr_r(r,sim)            cagr over the enitre time frame - total sectoral production for the domestic market - region r under physical shock sim
        yt_tcagr_us(sim)             cagr over the enitre time frame - total sectoral production for the domestic market - time period tp under physical shock sim

        yd_tcagr_sa(r,s_a,sim)       cagr over the enitre time frame - sectoral production for the domestic market - region r for aggregate sector s_a at time period tp under physical shock sim
        yxd_tcagr_sa(r,s_a,sim)      cagr over the enitre time frame - sectoral production for the intranational market - region r for aggregate sector s_a at time period tp under physical shock sim
        yxf_tcagr_sa(r,s_a,sim)      cagr over the enitre time frame - sectoral production for the international market - region r for aggregate sector s_a at time period tp under physical shock sim
        yt_tcagr_sa(r,s_a,sim)       cagr over the enitre time frame - total sectoral production for the domestic market - region r for aggregate sector s_a at time period tp under physical shock sim

        grpfc_tcagr_r(r,sim)         cagr over the enitre time frame - gross regional product (final demand) - region r under physical shock sim
        grpfc_tcagr_us(sim)          cagr over the enitre time frame - gross domestic product (final demand) - physical shock sim
        grpva_tcagr_r(r,sim)         cagr over the enitre time frame - gross regional product (value added) - region r under physical shock sim
        grpva_tcagr_us(sim)          cagr over the enitre time frame - gross domestic poduct (value added - physical shock sim
;

* ---------------------------------------------------------------------------- *
*        Load population step 2010 - 2050 (in 5 year intervals):
* ---------------------------------------------------------------------------- *

$GDXIN %load_pop%
$LOAD us_pop_proj
$GDXIN

$GDXIN %load_pop_carg%
$LOAD us_pop_proj_cagr
$GDXIN

poprate(r,simtp) =  us_pop_proj_cagr(r,simtp);

        DISPLAY   us_pop_proj, us_pop_proj_cagr, poprate;

* ---------------------------------------------------------------------------- *
*        Load macroeconomic variables:
* ---------------------------------------------------------------------------- *

*               Load macroeconomic parameters from sensitivity analysis

$GDXIN %load_macro%
$LOAD drate irate kprate lprate
$GDXIN

        DISPLAY drate, irate, kprate, lprate;

* ---------------------------------------------------------------------------- *
*        WDRG coupling:
* ---------------------------------------------------------------------------- *

*               Label for Demand Response coupling

$label wdrg_coupling

$if not %s_wd%==1 $goto demand_response_coupling

*               Run WDRG pre-process code

$include %load_wdrg_pre%

PARAMETER
        wdrgtfp_adj(r,crps,sim,tp)  TFP adjustment parameters for the agricultural yield impacts
        wdrgelas_adj                Parameters used to change the elasticities of transformation and trade to 0 for the trouble crop sectors
;

*       Load WDRG Yield Impacts

*$GDXIN %load_yield_ot%
*$load inverse_yld
*$GDXIN

*               Re-define WDRG TFP

wdrgtfp_adj(r,crps,forloop,simtp) = 1.00;
wdrgtfp_adj(r,crps,'bau','2015')  = inverse_yld(crps,r);
*wdrgtfp_adj('ROUS',crps,'bau',simtp)      = 1.00;
*wdrgtfp_adj('ROUS',crps,'cntf_wdrg',simtp) = 1.00;

*               Re-define elasticities for crop sectors

wdrgelas_adj('eta',crps)  =  1.00;
wdrgelas_adj('top',crps)  =  1.00;
wdrgelas_adj('bot',crps)  =  1.00;

wdrgelas_adj('eta','cba') =  0.25;
wdrgelas_adj('top','cba') =  0.25;
wdrgelas_adj('bot','cba') =  0.10;

wdrgelas_adj('eta','grn') =  0.25;
wdrgelas_adj('top','grn') =  0.25;
wdrgelas_adj('bot','grn') =  0.10;

wdrgelas_adj('eta','pfb') =  0.25;
wdrgelas_adj('top','pfb') =  0.25;
wdrgelas_adj('bot','pfb') =  0.10;

wdrgelas_adj('eta','osa') =  0.00;
wdrgelas_adj('top','osa') =  0.00;
wdrgelas_adj('bot','osa') =  0.00;

wdrgelas_adj('eta','oca') =  0.00;
wdrgelas_adj('top','oca') =  0.00;
wdrgelas_adj('bot','oca') =  0.00;

et_dx(crps) = et_dx(crps) * wdrgelas_adj('eta',crps);
es_nf(crps) = es_nf(crps) * wdrgelas_adj('top',crps);
es_dn(crps) = es_dn(crps) * wdrgelas_adj('bot',crps);

    Display wdrgtfp_adj, wdrgelas_adj;

* ---------------------------------------------------------------------------- *
*        Demand Response coupling:
* ---------------------------------------------------------------------------- *

*               Label for Demand Response coupling

$label demand_response_coupling

$if not %s_dr%==1 $goto psm_nonproductivity_coupling

*               Load demand response GDX from R

PARAMETER
        pchg_sts(r)                         percentage change in per capita load from the baseline - demand response from space-for-time substitution - region r
        pchg_td(r)                          percentage change in per capita load from the baseline - demand response from temporal downscaling - region r
        shk_dr_pchg_prod(r,g,s,tp,sim)      percentage change in per capita load from the baseline - demand response - region r for sector s at time period tp under physical shock sim
        shk_dr_pchg_consum(r,s,h,tp,sim)    percentage change in per capita load from the baseline - demand response - region r for household h at time period tp under physical shock sim
        dre_prod_sim(r,g,s,tp,sim)          Electricity Demand response productivity shock by sector and time period
        dre_hh_sim(r,s,h,tp,sim)            Electricity Demand response productivity shock by household and time period
;

$GDXIN %load_drm%
* $LOAD pchg_sts pchg_td
$LOAD pchg_sts
$GDXIN

*               Define demand response shock

* space time substitution
shk_dr_pchg_prod(r,"ele",s,tlast,sim_drm_sts)   = pchg_sts(r);
shk_dr_pchg_consum(r,"ele",h,tlast,sim_drm_sts) = pchg_sts(r);

* temporal downscaling -- not using
*shk_dr_pchg_prod(r,"ele",s,tlast,sim_drm_td)    = pchg_td(r);
*shk_dr_pchg_consum(r,"ele",h,tlast,sim_drm_td)  = pchg_td(r);

*        DISPLAY pchg_sts, pchg_td, shk_dr_pchg_prod, shk_dr_pchg_consum;
        DISPLAY pchg_sts, shk_dr_pchg_prod, shk_dr_pchg_consum;

* ---------------------------------------------------------------------------- *
*        Power System Model without productivity impacts coupling:
* ---------------------------------------------------------------------------- *

*               Label for PSM coupling (no productivity impacts)

$label psm_nonproductivity_coupling

$if not %s_pw%==1 $goto psm_productivity_coupling

*               Load PSM TFP scalars from PSM:

PARAMETER
        ele_impacts(*)            PSM electricity productivity impact for economic region r
        psm_ele_pi(r)
        psmtfp_temp(r,g,s,tp,sim)     Temporary PSM productivity impact for good g in region r
        psmtfp_fa_temp(r,s,tp,sim)    Temporary PSM productivity impact for sector s in region r
;

$GDXIN %load_psm_ele%
$LOAD ele_impacts
$GDXIN

*               Define PSM Shock for electricity in each region
loop(econr_map(econr,r),
         psm_ele_pi(r) =  ele_impacts(econr);
);

*amc_chg('CA')      = ChangeCali;
*amc_chg('ROWECC')  = ChangeRO;

*amc_chg('CA')      = 0.029944406003197;  -- PSM with DR first iteration
*amc_chg('ROWECC')  = 0.0129947901071921; -- PSM with DR first iteration

*ws
*amc_chg('CA')      = 0.001992006;
*amc_chg('ROWECC')  = 0.008566497;

*ret
*amc_chg('CA')      = 0.084982763;
*amc_chg('ROWECC')  = 0.079356422;

*ws+ret
*amc_chg('CA') = 0.079371377;
*amc_chg('ROWECC') = 0.083523413;

*psm_epi('cal') = 0.079371377;
*psm_epi('arz') = 0.083523413;
*psm_epi('col') = 0.083523413;
*psm_epi('ido') = 0.083523413;
*psm_epi('mta') = 0.083523413;
*psm_epi('nmo') = 0.083523413;
*psm_epi('nva') = 0.083523413;
*psm_epi('org') = 0.083523413;
*psm_epi('uth') = 0.083523413;
*psm_epi('was') = 0.083523413;
*psm_epi('wyo') = 0.083523413;

psmtfp_temp(r,g,s,simtp,forloop)  = 1;
psmtfp_fa_temp(r,s,simtp,forloop) = 1;

psmtfp_temp(r,g,'ele',tlast,sim_psm)  = 1 + psm_ele_pi(r);
psmtfp_fa_temp(r,'ele',tlast,sim_psm) = 1 + psm_ele_pi(r);

*        DISPLAY ChangeCali, ChangeRO, amc_chg, psmtfp_fa_temp, psmtfp_temp;
        DISPLAY psm_ele_pi, psmtfp_fa_temp, psmtfp_temp;

* ---------------------------------------------------------------------------- *
*        Power System Model with productivity impacts coupling:
* ---------------------------------------------------------------------------- *

*               Label for PSM coupling (no productivity impacts)

$label psm_productivity_coupling

$if not %s_pf%==1 $goto psm_nosul_multi_coupling

*               Load Sullivan Impacts
PARAMETER
    sull_impacts(*,*)     Change in production cost from outages (see "Sullivan electricity loss calculations.xls")
    psm_nele_pi(s,r)
;

$GDXIN %load_psm_nele%
$load sull_impacts
$GDXIN

*               Define PSM production cost (Sullivan production loss calculations to PSM TFP)
loop(econr_map(econr,r),
         psm_nele_pi(nele,r) = sull_impacts(nele,econr);
);
* psmtfp_temp(cal,g,nele,simtp,sim_psmfull)  = 1 + Sull_impacts(nele,'CA');
* psmtfp_fa_temp(cal,nele,simtp,sim_psmfull) = 1 + Sull_impacts(nele,'CA');

* psmtfp_temp(rwc,g,nele,tlast,sim_psmfull)  = 1 + Sull_impacts(nele,'ROWECC');
* psmtfp_fa_temp(rwc,nele,tlast,sim_psmfull) = 1 + Sull_impacts(nele,'ROWECC');

psmtfp_temp(r,g,nele,simtp,sim_psmfull)  = 1 + psm_nele_pi(nele,r);
psmtfp_fa_temp(r,nele,simtp,sim_psmfull) = 1 + psm_nele_pi(nele,r);

        DISPLAY sull_impacts, psm_nele_pi, psmtfp_fa_temp, psmtfp_temp;

*               Include the PSM electicity simulation in when the productiviey impacts are turned on

*psmtfp_temp(cal,g,'ele',tlast,sim_psmfull)  = 1 + amc_chg('CA');
*psmtfp_fa_temp(cal,'ele',tlast,sim_psmfull) = 1 + amc_chg('CA');

*psmtfp_temp(rwc,g,'ele',tlast,sim_psmfull)  = 1 + amc_chg('ROWECC');
*psmtfp_fa_temp(rwc,'ele',tlast,sim_psmfull) = 1 + amc_chg('ROWECC');

psmtfp_temp(r,g,'ele',tlast,sim_psmfull)  = 1 + psm_ele_pi(r);
psmtfp_fa_temp(r,'ele',tlast,sim_psmfull) = 1 + psm_ele_pi(r);

*        DISPLAY ChangeCali, ChangeRO, amc_chg, psmtfp_fa_temp, psmtfp_temp;
        DISPLAY psm_ele_pi, psmtfp_fa_temp, psmtfp_temp;

* ---------------------------------------------------------------------------- *
*        Power System Model with other coupled models:
* ---------------------------------------------------------------------------- *

$label psm_nosul_multi_coupling

$if not %s_pw%==1 $goto psm_full_multi_coupling

*               Include the PSM electicity simulation in when the productiviey impacts are turned on

*psmtfp_temp(cal,g,'ele',tlast,policy)  = 1 + amc_chg('CA');
*psmtfp_fa_temp(cal,'ele',tlast,policy) = 1 + amc_chg('CA');

*psmtfp_temp(rwc,g,'ele',tlast,policy)  = 1 + amc_chg('ROWECC');
*psmtfp_fa_temp(rwc,'ele',tlast,policy) = 1 + amc_chg('ROWECC');

psmtfp_temp(r,g,'ele',tlast,policy)  = 1 + psm_ele_pi(r);
psmtfp_fa_temp(r,'ele',tlast,policy) = 1 + psm_ele_pi(r);

        DISPLAY psm_ele_pi, psmtfp_fa_temp, psmtfp_temp;

$label psm_full_multi_coupling

* ---------------------------------------------------------------------------- *
*        Power System Model productivity impacts with other coupled models:
* ---------------------------------------------------------------------------- *

$if not %s_pf%==1 $goto psm_dr_coupling

*               Define PSM production cost (Sullivan production loss calculations to PSM TFP)

*psmtfp_temp(cal,g,nele,tlast,policy)  = 1 + Sull_impacts(nele,'CA');
*psmtfp_fa_temp(cal,nele,tlast,policy) = 1 + Sull_impacts(nele,'CA');

*psmtfp_temp(rwc,g,nele,tlast,policy)  = 1 + Sull_impacts(nele,'ROWECC');
*psmtfp_fa_temp(rwc,nele,tlast,policy) = 1 + Sull_impacts(nele,'ROWECC');

psmtfp_temp(r,g,nele,tlast,policy)  = 1 + psm_nele_pi(nele,r);
psmtfp_fa_temp(r,nele,tlast,policy) = 1 + psm_nele_pi(nele,r);

        DISPLAY psm_nele_pi, psmtfp_fa_temp, psmtfp_temp;

* ---------------------------------------------------------------------------- *
*        Demand Response Emulator impacts with other coupled models:
* ---------------------------------------------------------------------------- *

$label psm_dr_coupling

$if not %s_dr%==1 $goto recursive_loop

*               Define demand response shock

* space time substitution
shk_dr_pchg_prod(r,"ele",s,tlast,policy)   = pchg_sts(r);
shk_dr_pchg_consum(r,"ele",h,tlast,policy) = pchg_sts(r);

* temporal downscaling -- not using
*shk_dr_pchg_prod(r,"ele",s,tlast,policy)    = pchg_td(r);
*shk_dr_pchg_consum(r,"ele",h,tlast,policy)  = pchg_td(r);

*        DISPLAY pchg_sts, pchg_td, shk_dr_pchg_prod, shk_dr_pchg_consum;
        DISPLAY pchg_sts;

* ---------------------------------------------------------------------------- *
*        Begin Policy Simulation Recursive Dynamic Loop:
* ---------------------------------------------------------------------------- *

*               Label for recursive loop

$label recursive_loop

loop(forloop(sim),

* -------------------------------------------- *
* Define intial conditions
* -------------------------------------------- *
        tfp(r,s)          = 1.00;
        dre_prod(r,g,s)   = 1.00;
        dre_hh(r,s,h)     = 1.00;
        psmtfp(r,g,s)     = 1.00;
        psmtfp_fa(r,s)    = 1.00;
        lpidx(r,tp)       = 1.00;
        vdmi(r,s)         = vdmi0(r,s);
        vdgm(r,s,pub)     = vdgm0(r,s,pub);
        vigm(r,s,trd,pub) = vigm0(r,s,trd,pub);
        vipm(r,s,trd,h)   = vipm0(r,s,trd,h);
        vdifm(r,g,s)      = vdifm0(r,g,s);
        vfm(r,fa,s)       = vfm0(r,fa,s);
        evo(r,i,f)        = evo0(r,i,f);
        vfm_k(r,s)        = vfm_k0(r,s);
        vfm_l(r,s)        = vfm_l0(r,s);
        evo_k(r,h)        = evo_k0(r,h);
        evo_l(r,h)        = evo_l0(r,h);
        va(r,s)           = va0(r,s);
        vn(g)             = vn0(g) ;
        vinv(r)           = vinv0(r);
        vinvd(r,g)        = vinvd0(r,g);
        vinvh(r,h)        = vinvh0(r,h);
        incadj(r,h)       = incadj0(r,h);
        vxm(r,g,trd)      = vxm0(r,g,trd);

*   Define base year wage payments (or capital payments) and labor endowments for household
        evoi_k0(r,h)   = evo_k0(r,h);
        evoi_l0(r,h)   = evo_l0(r,h);

*   Define base year capital and labor stock
        kstock0(r)    = sum(h, evoi_k0(r,h))/(drate+irate);
        lstock0(r)    = sum(h, evoi_l0(r,h));

*   Define base year population
*        us_pop_proj0(r) = us_pop_proj(r,'2010');

*   Define inital time period capital and labor stock
        kstock(r,tfirst) = kstock0(r);
        lstock(r,tfirst) = lstock0(r);

* -------------------------------------------- *
*   Begin Recursive loop
* -------------------------------------------- *
        loop(simtp,

*               Define the WDRG total factor productivities for crop sectors
$label define_wdrg
$if not %s_wd%==1 $goto define_drm

                wdrgtfp(r,crps) = wdrgtfp_adj(r,crps,sim,simtp);

*               Define the demand response shock for the production and utility functions
$label define_drm
$if not %s_dr%==1 $goto define_psmfp

                 dre_prod(r,"ele",s) = 1 + shk_dr_pchg_prod(r,'ele',s,simtp,sim);
                 dre_hh(r,"ele",h)   = 1 + shk_dr_pchg_consum(r,'ele',h,simtp,sim);

*                Check that the percentage change shock is being read into the model correctly
                 shk_dr_pchg_prod_sim(r,g,s,simtp,sim)   = shk_dr_pchg_prod(r,g,s,simtp,sim);
                 shk_dr_pchg_consum_sim(r,s,h,simtp,sim) = shk_dr_pchg_consum(r,s,h,simtp,sim);

*               Check that the demand response shock is being calculated correctly
                dre_prod_sim(r,g,s,simtp,sim) = dre_prod(r,g,s);
                dre_hh_sim(r,s,h,simtp,sim)   = dre_hh(r,s,h);

*               Define PSM productivity impacts on electricity sector and
*                      PSM production cost (Sullivan impacts) on non-electricity sectors (if defined above)
$label define_psmfp
$if not %s_pw%==1 $goto call_mpsge

                psmtfp(r,g,s)  = psmtfp_temp(r,g,s,simtp,sim);
                psmtfp_fa(r,s) = psmtfp_fa_temp(r,s,simtp,sim);

*               Save PSM productivity shocks
                psmtfp_sim(r,g,s,simtp,sim)  =  psmtfp(r,g,s);
                psmtfp_fa_sim(r,s,simtp,sim) =  psmtfp_fa(r,s);

*               Call MPSGE solve
$label call_mpsge

                soe.iterlim = 1000000;
                calib = 0;
$include SOE.GEN
                solve soe using mcp;

*       Define the recusive parameter
$if %t_y%==2010 $goto model_parameters

*               Define population index in the next time period
                lpidx(r,simtp+1) = lpidx(r,simtp)*(1 + poprate(r,simtp+1) + lprate)**ts(simtp+1);

*               Define employment stock in the next time period
                lstock(r,simtp+1)= lstock0(r)*lpidx(r,simtp+1);

*               Define capital stock (capital acculation) in the next time period
*                kstock(r,simtp+1)= ts(simtp+1) * sum(hh, rhscale.l(r,hh)*vinvh(r,hh)) + (1 - drate)**ts(simtp+1) * kstock(r,simtp);
                kstock(r,simtp+1)= ts(simtp+1) * vinv(r)*inv.l(r) + (1 - drate)**ts(simtp+1) * kstock(r,simtp);

*               Capital stock productivity
                kstock(r,simtp+1)= kstock(r,simtp+1)*(1 + kprate);

*               Translate employment and capital stocks back to MPSGE Model
                evo_k(r,h) = evoi_k0(r,h)*kstock(r,simtp+1)/kstock0(r);
                evo_l(r,h) = evoi_l0(r,h)*lstock(r,simtp+1)/lstock0(r);

                incadj_sim(r,h,simtp,sim)  = incadj(r,h);
                rhscale_sim(r,h,simtp,sim) = rhscale.l(r,h);
                evolp_sim(r,h,simtp,sim)   = evolp.l(r,h);
                evokp_sim(r,h,simtp,sim)   = evokp.l(r,h);
                ca_sim(r,s,h,tp,sim)       = ca.l(r,s,h);
                hinv_sim(r,h,tp,sim)       = hinv.l(r,h);
                incadj_cal_sim(r,h,tp,sim) = sum(s,ca.l(r,s,h)) + hinv.l(r,h)
                                                 + evokp.l(r,h) + evolp.l(r,h);

*        Save model output for each time period and policy forloop

*$label save_recursive_parameterss

*       Recursive stock and flows
*       Population Index
        lpidx_sim(r,simtp)          = lpidx(r,simtp);

*       Capital and Labor stock
        kstock_sim(r,simtp,sim)     = kstock(r,simtp);
        lstock_sim(r,simtp,sim)     = lstock(r,simtp);

*       Capital and Labor Flow
        evok_sim(r,h,simtp,sim)     = evoi_k0(r,h)*kstock(r,simtp+1)/kstock0(r);
        evol_sim(r,h,simtp,sim)     = evoi_l0(r,h)*lstock(r,simtp+1)/lstock0(r);

$label model_parameters

*       Price
        py_sim(r,s,simtp,sim)       = p.l(r,s);
        pa_sim(r,s,simtp,sim)       = pa.l(r,s);
        pc_sim(r,h,simtp,sim)       = pc.l(r,h);
        pn_sim(s,simtp,sim)         = pn.l(s);
        pinv_sim(r,simtp,sim)       = pinv.l(r);
        pgov_sim(r,pub,simtp,sim)   = pgov.l(r,pub);
        pk_sim(r,simtp,sim)         = pf_k.l(r);
        pl_sim(r,simtp,sim)         = pf_l.l(r);
        pbt_sim(r,simtp,sim)        = ptax.l(r);

*       Activity
        yd_sim(r,s,simtp,sim)       = yd.l(r,s);
        yd_sim_r(r,simtp,sim)       = sum(s, yd.l(r,s));
        yd_sim_us(simtp,sim)        = sum((r,s), yd.l(r,s));
        vyd_sim(r,s,simtp,sim)      = p.l(r,s)*yd.l(r,s);
        yxd_sim(r,s,simtp,sim)      = yxd.l(r,s);
        yxd_sim_r(r,simtp,sim)      = sum(s, yxd.l(r,s));
        yxd_sim_us(simtp,sim)       = sum((r,s), yxd.l(r,s));
        yxf_sim(r,s,simtp,sim)      = yxf.l(r,s);
        yxf_sim_r(r,simtp,sim)      = sum(s, yxf.l(r,s));
        yxf_sim_us(simtp,sim)       = sum((r,s), yxf.l(r,s));
        yt_sim(r,s,simtp,sim)       = yd.l(r,s) + yxd.l(r,s) + yxf.l(r,s);
        vyt_sim(r,s,simtp,sim)      = p.l(r,s)*yd.l(r,s) + pn.l(s)*yxd.l(r,s) + pfx.l*yxf.l(r,s);
        yt_sim_r(r,simtp,sim)       = sum(s, yd.l(r,s) + yxd.l(r,s) + yxf.l(r,s));
        yt_sim_us(simtp,sim)        = sum((r,s), yd.l(r,s) + yxd.l(r,s) + yxf.l(r,s));
        a_sim(r,s,simtp,sim)        = aa.l(r,s);

*       Activity by subset
        yd_sim_sa(r,s_a,simtp,sim)  = sum(s_map(s_a,s), yd.l(r,s));
        yxd_sim_sa(r,s_a,simtp,sim) = sum(s_map(s_a,s), yxd.l(r,s));
        yxf_sim_sa(r,s_a,simtp,sim) = sum(s_map(s_a,s), yxf.l(r,s));
        yt_sim_sa(r,s_a,simtp,sim)  = sum(s_map(s_a,s), yd.l(r,s) + yxd.l(r,s) + yxf.l(r,s));

        yt_share_sa(r,s,simtp,sim)$sum(s_map(s_a,s), yd.l(r,s) + yxd.l(r,s) + yxf.l(r,s)) = (yd.l(r,s) + yxd.l(r,s) + yxf.l(r,s))/sum(s_map(s_a,s), yd.l(r,s) + yxd.l(r,s) + yxf.l(r,s));
        yt_share_r(r,s,simtp,sim)$(sum(rr, yd.l(r,s) + yxd.l(r,s) + yxf.l(r,s)))          = (yd.l(r,s) + yxd.l(r,s) + yxf.l(r,s))/sum(rr, yd.l(r,s) + yxd.l(r,s) + yxf.l(r,s));

*       Demand
        kd_sim(r,s,simtp,sim)       = kd.l(r,s);
        kd_sim_r(r,simtp,sim)       = sum(s, kd.l(r,s));
        kd_sim_us(simtp,sim)        = sum((r,s), kd.l(r,s));
        ld_sim(r,s,simtp,sim)       = ld.l(r,s);
        ld_sim_r(r,simtp,sim)       = sum(s, ld.l(r,s));
        ld_sim_us(simtp,sim)        = sum((r,s), ld.l(r,s));
        btd_sim(r,s,simtp,sim)      = btd.l(r,"btax",s);
        btd_sim_r(r,simtp,sim)      = sum(s, btd.l(r,"btax",s));
        btd_sim_us(simtp,sim)       = sum((r,s), btd.l(r,"btax",s));
        id_sim(r,g,s,simtp,sim)     = id.l(r,g,s);
        id_sim_g(g,simtp,sim)       = sum((r,s), id.l(r,g,s));
        id_sim_r(r,simtp,sim)       = sum((g,s), id.l(r,g,s));
        id_sim_us(simtp,sim)        = sum((r,g,s), id.l(r,g,s));

        md_sim(r,s,simtp,sim)       = md.l(r,s);
        md_sim_r(r,simtp,sim)       = sum(s, md.l(r,s));
        md_sim_us(simtp,sim)        = sum((r,s), md.l(r,s));
        mf_sim(r,s,simtp,sim)       = mf.l(r,s);
        mf_sim_r(r,simtp,sim)       = sum(s, mf.l(r,s));
        mf_sim_us(simtp,sim)        = sum((r,s), mf.l(r,s));

        cd_sim(r,h,simtp,sim)       = cd.l(r,h);
        cd_sim_h(h,simtp,sim)       = sum(r, cd.l(r,h));
        cd_sim_r(r,simtp,sim)       = sum(h, cd.l(r,h));
        cd_sim_us(simtp,sim)        = sum((r,h), cd.l(r,h));

        i_sim_r(r,simtp,sim)        = invd.l(r);
        i_sim_us(simtp,sim)         = sum(r, invd.l(r));

        g_sim(r,pub,simtp,sim)      = govd.l(r,pub);
        g_sim_r(r,simtp,sim)        = sum(pub, govd.l(r,pub));
        g_sim_us(simtp,sim)         = sum((r,pub), govd.l(r,pub));

        trev_sim_r(r,simtp,sim)     = taxrev.l(r);
        trev_sim_us(simtp,sim)      = sum(r, taxrev.l(r));

        grpfc_sim_r(r,simtp,sim)    = sum(h, cd.l(r,h)) + sum(pub, govd.l(r,pub)) + invd.l(r) - sum(s, mf.l(r,s)) +  sum(s, yxf.l(r,s));
        grpfc_sim_us(simtp,sim)     = sum((r,h), cd.l(r,h)) + sum((r,pub), govd.l(r,pub)) + sum(r, invd.l(r)) - sum((r,s), mf.l(r,s)) +  sum((r,s), yxf.l(r,s));
        grpva_sim_r(r,simtp,sim)    = pf_k.l(r)*sum(s, kd.l(r,s)) + pf_l.l(r)*sum(s, ld.l(r,s)) + ptax.l(r)*sum(s, btd.l(r,"btax",s));
        grpva_sim_us(simtp,sim)     = sum(r, pf_k.l(r)*sum(s, kd.l(r,s))) + sum(r, pf_l.l(r)*sum(s, ld.l(r,s))) + sum(r, ptax.l(r)*sum(s, btd.l(r,"btax",s)));
    );
);

* ---------------------------------------------------------------------------- *
*        Between period compound annualized growth rate:
* ---------------------------------------------------------------------------- *

$label year_growth_rate
$if %t_y%==2010 $goto total_growth_rate

loop(forloop(sim),
        loop(simtp,

*           Activity
            yd_cagr(r,s,simtp+1,sim)$(yd_sim(r,s,simtp,sim)<> 0)   = (yd_sim(r,s,simtp+1,sim)/yd_sim(r,s,simtp,sim))**(1/ts(simtp+1))-1;
            yd_cagr_r(r,simtp+1,sim)$(yd_sim_r(r,simtp,sim)<> 0)   = (yd_sim_r(r,simtp+1,sim)/yd_sim_r(r,simtp,sim))**(1/ts(simtp+1))-1;
            yd_cagr_us(simtp+1,sim)$(yd_sim_us(simtp,sim)<> 0)     = (yd_sim_us(simtp+1,sim)/yd_sim_us(simtp,sim))**(1/ts(simtp+1))-1;
            yxd_cagr(r,s,simtp+1,sim)$(yxd_sim(r,s,simtp,sim)<> 0) = (yxd_sim(r,s,simtp+1,sim)/yxd_sim(r,s,simtp,sim))**(1/ts(simtp+1))-1;
            yxd_cagr_r(r,simtp+1,sim)$(yxd_sim_r(r,simtp,sim)<> 0) = (yxd_sim_r(r,simtp+1,sim)/yxd_sim_r(r,simtp,sim))**(1/ts(simtp+1))-1;
            yxd_cagr_us(simtp+1,sim)$(yxd_sim_us(simtp,sim)<> 0)   = (yxd_sim_us(simtp+1,sim)/yxd_sim_us(simtp,sim))**(1/ts(simtp+1))-1;
            yxf_cagr(r,s,simtp+1,sim)$(yxf_sim(r,s,simtp,sim)<> 0) = (yxf_sim(r,s,simtp+1,sim)/yxf_sim(r,s,simtp,sim))**(1/ts(simtp+1))-1;
            yxf_cagr_r(r,simtp+1,sim)$(yxf_sim_r(r,simtp,sim)<> 0) = (yxf_sim_r(r,simtp+1,sim)/yxf_sim_r(r,simtp,sim))**(1/ts(simtp+1))-1;
            yxf_cagr_us(simtp+1,sim)$(yxf_sim_us(simtp,sim)<> 0)   = (yxf_sim_us(simtp+1,sim)/yxf_sim_us(simtp,sim))**(1/ts(simtp+1))-1;
            yt_cagr(r,s,simtp+1,sim)$(yt_sim(r,s,simtp,sim)<> 0)   = (yt_sim(r,s,simtp+1,sim)/yt_sim(r,s,simtp,sim))**(1/ts(simtp+1))-1;
            yt_cagr_r(r,simtp+1,sim)$(yt_sim_r(r,simtp,sim)<> 0)   = (yt_sim_r(r,simtp+1,sim)/yt_sim_r(r,simtp,sim))**(1/ts(simtp+1))-1;
            yt_cagr_us(simtp+1,sim)$(yt_sim_us(simtp,sim)<> 0)     = (yt_sim_us(simtp+1,sim)/yt_sim_us(simtp,sim))**(1/ts(simtp+1))-1;

*           Activity by aggregate sector
            yd_cagr_sa(r,s_a,simtp+1,sim)$(yd_sim_sa(r,s_a,simtp,sim)<> 0)   = (yd_sim_sa(r,s_a,simtp+1,sim)/yd_sim_sa(r,s_a,simtp,sim))**(1/ts(simtp+1))-1;
            yxd_cagr_sa(r,s_a,simtp+1,sim)$(yxd_sim_sa(r,s_a,simtp,sim)<> 0) = (yxd_sim_sa(r,s_a,simtp+1,sim)/yxd_sim_sa(r,s_a,simtp,sim))**(1/ts(simtp+1))-1;
            yxf_cagr_sa(r,s_a,simtp+1,sim)$(yxf_sim_sa(r,s_a,simtp,sim)<> 0) = (yxf_sim_sa(r,s_a,simtp+1,sim)/yxf_sim_sa(r,s_a,simtp,sim))**(1/ts(simtp+1))-1;
            yt_cagr_sa(r,s_a,simtp+1,sim)$(yt_sim_sa(r,s_a,simtp,sim)<> 0)   = (yt_sim_sa(r,s_a,simtp+1,sim)/yt_sim_sa(r,s_a,simtp,sim))**(1/ts(simtp+1))-1;

*           GRP
            grpva_cagr_r(r,simtp+1,sim) = (grpva_sim_r(r,simtp+1,sim)/grpva_sim_r(r,simtp,sim))**(1/ts(simtp+1))-1;
            grpva_cagr_us(simtp+1,sim)  = (grpva_sim_us(simtp+1,sim)/grpva_sim_us(simtp,sim))**(1/ts(simtp+1))-1;
            grpfc_cagr_r(r,simtp+1,sim) = (grpfc_sim_r(r,simtp+1,sim)/grpfc_sim_r(r,simtp,sim))**(1/ts(simtp+1))-1;
            grpfc_cagr_us(simtp+1,sim)  = (grpfc_sim_us(simtp+1,sim)/grpfc_sim_us(simtp,sim))**(1/ts(simtp+1))-1;

        );
);

* ---------------------------------------------------------------------------- *
*        Compound annualized growth rate other whole time frame:
* ---------------------------------------------------------------------------- *

$label total_growth_rate
$if %t_y%==2010 $goto pchg_calculation

loop(forloop(sim),

*       Activity
        yd_tcagr(r,s,sim)$(yd_sim(r,s,'2010',sim) <> 0)    = (yd_sim(r,s,'%t_y%',sim)/yd_sim(r,s,'2010',sim))**(1/(tts))-1;
        yd_tcagr_r(r,sim)                                  = (yd_sim_r(r,'%t_y%',sim)/yd_sim_r(r,'2010',sim))**(1/(tts))-1;
        yd_tcagr_us(sim)                                   = (yd_sim_us('%t_y%',sim)/yd_sim_us('2010',sim))**(1/(tts))-1;
        yxd_tcagr(r,s,sim)$(yxd_sim(r,s,'2010',sim) <> 0)  = (yxd_sim(r,s,'%t_y%',sim)/yxd_sim(r,s,'2010',sim))**(1/(tts))-1;
        yxd_tcagr_r(r,sim)                                 = (yxd_sim_r(r,'%t_y%',sim)/yxd_sim_r(r,'2010',sim))**(1/(tts))-1;
        yxd_tcagr_us(sim)                                  = (yxd_sim_us('%t_y%',sim)/yxd_sim_us('2010',sim))**(1/(tts))-1;
        yxf_tcagr(r,s,sim)$(yxf_sim(r,s,'2010',sim) <> 0)  = (yxf_sim(r,s,'%t_y%',sim)/yxf_sim(r,s,'2010',sim))**(1/(tts))-1;
        yxf_tcagr_r(r,sim)                                 = (yxf_sim_r(r,'%t_y%',sim)/yxf_sim_r(r,'2010',sim))**(1/(tts))-1;
        yxf_tcagr_us(sim)                                  = (yxf_sim_us('%t_y%',sim)/yxf_sim_us('2010',sim))**(1/(tts))-1;
        yt_tcagr(r,s,sim)$(yt_sim(r,s,'2010',sim) <> 0)    = (yt_sim(r,s,'%t_y%',sim)/yt_sim(r,s,'2010',sim))**(1/(tts))-1;
        yt_tcagr_r(r,sim)                                  = (yt_sim_r(r,'%t_y%',sim)/yt_sim_r(r,'2010',sim))**(1/(tts))-1;
        yt_tcagr_us(sim)                                   = (yt_sim_us('%t_y%',sim)/yt_sim_us('2010',sim))**(1/(tts))-1;

*       Activity by aggregate sector
        yd_tcagr_sa(r,s_a,sim)$(yd_sim_sa(r,s_a,'2010',sim) <> 0)   = (yd_sim_sa(r,s_a,'%t_y%',sim)/yd_sim_sa(r,s_a,'2010',sim))**(1/(tts))-1;
        yxd_tcagr_sa(r,s_a,sim)$(yxd_sim_sa(r,s_a,'2010',sim) <> 0) = (yxd_sim_sa(r,s_a,'%t_y%',sim)/yxd_sim_sa(r,s_a,'2010',sim))**(1/(tts))-1;
        yxf_tcagr_sa(r,s_a,sim)$(yxf_sim_sa(r,s_a,'2010',sim) <> 0) = (yxf_sim_sa(r,s_a,'%t_y%',sim)/yxf_sim_sa(r,s_a,'2010',sim))**(1/(tts))-1;
        yt_tcagr_sa(r,s_a,sim)$(yt_sim_sa(r,s_a,'2010',sim) <> 0)   = (yt_sim_sa(r,s_a,'%t_y%',sim)/yt_sim_sa(r,s_a,'2010',sim))**(1/(tts))-1;

*       GRP
        grpfc_tcagr_r(r,sim) = (grpfc_sim_r(r,'%t_y%',sim)/grpfc_sim_r(r,'2010',sim))**(1/(tts))-1;
        grpfc_tcagr_us(sim)  = (grpfc_sim_us('%t_y%',sim)/grpfc_sim_us('2010',sim))**(1/(tts))-1;
        grpva_tcagr_r(r,sim) = (grpva_sim_r(r,'%t_y%',sim)/grpva_sim_r(r,'2010',sim))**(1/(tts))-1;
        grpva_tcagr_us(sim)  = (grpva_sim_us('%t_y%',sim)/grpva_sim_us('2010',sim))**(1/(tts))-1;

);

* ---------------------------------------------------------------------------- *
*        Percentage Change from baseline (for period 9 only):
* ---------------------------------------------------------------------------- *

$label pchg_calculation

*               Define policy loop for precentage change

SET
* Dynamic subsets
      policy_pchg(sim)   Simulations for calculating percentage change
;

* include the policy and policy BAU simulations
policy_pchg(sim)$(policy(sim))     = YES;
policy_pchg(sim)$(policy_bau(sim)) = YES;

* exclude the BAU simulation
policy_pchg(sim)$(sim_bau(sim))    = NO;

        DISPLAY  policy_pchg;

LOOP(policy(sim),
    LOOP(policy_bau(ssim),
*       Price
        py_pchg(r,s,policy_pchg,policy_bau)$(py_sim(r,s,'%t_y%',policy_bau) <> 0)         = (py_sim(r,s,'%t_y%',policy_pchg)/py_sim(r,s,'%t_y%',policy_bau))-1;
        pa_pchg(r,s,policy_pchg,policy_bau)$(pa_sim(r,s,'%t_y%',policy_bau) <> 0)         = (pa_sim(r,s,'%t_y%',policy_pchg)/pa_sim(r,s,'%t_y%',policy_bau))-1;
        pc_pchg(r,h,policy_pchg,policy_bau)$(pc_sim(r,h,'%t_y%',policy_bau) <> 0)         = (pc_sim(r,h,'%t_y%',policy_pchg)/pc_sim(r,h,'%t_y%',policy_bau))-1;
        pn_pchg(s,policy_pchg,policy_bau)$(pn_sim(s,'%t_y%',policy_bau) <> 0)             = (pn_sim(s,'%t_y%',policy_pchg)/pn_sim(s,'%t_y%',policy_bau))-1;
        pinv_pchg(r,policy_pchg,policy_bau)$(pinv_sim(r,'%t_y%',policy_bau) <> 0)         = (pinv_sim(r,'%t_y%',policy_pchg)/pinv_sim(r,'%t_y%',policy_bau))-1;
        pgov_pchg(r,pub,policy_pchg,policy_bau)$(pgov_sim(r,pub,'%t_y%',policy_bau) <> 0) = (pgov_sim(r,pub,'%t_y%',policy_pchg)/pgov_sim(r,pub,'%t_y%',policy_bau))-1;
        pk_pchg(r,policy_pchg,policy_bau)$(pk_sim(r,'%t_y%',policy_bau) <> 0)             = (pk_sim(r,'%t_y%',policy_pchg)/pk_sim(r,'%t_y%',policy_bau))-1;
        pl_pchg(r,policy_pchg,policy_bau)$(pl_sim(r,'%t_y%',policy_bau) <> 0)             = (pl_sim(r,'%t_y%',policy_pchg)/pl_sim(r,'%t_y%',policy_bau))-1;
        pbt_pchg(r,policy_pchg,policy_bau)$(pbt_sim(r,'%t_y%',policy_bau) <> 0)           = (pbt_sim(r,'%t_y%',policy_pchg)/pbt_sim(r,'%t_y%',policy_bau))-1;

*       Activity
        yd_pchg(r,s,policy_pchg,policy_bau)$(yd_sim(r,s,'%t_y%',policy_bau) <> 0)         = (yd_sim(r,s,'%t_y%',policy_pchg)/yd_sim(r,s,'%t_y%',policy_bau))-1;
        yd_pchg_r(r,policy_pchg,policy_bau)$(yd_sim_r(r,'%t_y%',policy_bau) <> 0)         = (yd_sim_r(r,'%t_y%',policy_pchg)/yd_sim_r(r,'%t_y%',policy_bau))-1;
        yd_pchg_us(policy_pchg,policy_bau)$(yd_sim_us('%t_y%',policy_bau) <> 0)           = (yd_sim_us('%t_y%',policy_pchg)/yd_sim_us('%t_y%',policy_bau))-1;
        vyd_pchg(r,s,policy_pchg,policy_bau)$(vyd_sim(r,s,'%t_y%',policy_bau) <> 0)       = (vyd_sim(r,s,'%t_y%',policy_pchg)/vyd_sim(r,s,'%t_y%',policy_bau))-1;
        yxd_pchg(r,s,policy_pchg,policy_bau)$(yxd_sim(r,s,'%t_y%',policy_bau) <> 0)       = (yxd_sim(r,s,'%t_y%',policy_pchg)/yxd_sim(r,s,'%t_y%',policy_bau))-1;
        yxd_pchg_r(r,policy_pchg,policy_bau)$(yxd_sim_r(r,'%t_y%',policy_bau) <> 0)       = (yxd_sim_r(r,'%t_y%',policy_pchg)/yxd_sim_r(r,'%t_y%',policy_bau))-1;
        yxd_pchg_us(policy_pchg,policy_bau)$(yxd_sim_us('%t_y%',policy_bau) <> 0)         = (yxd_sim_us('%t_y%',policy_pchg)/yxd_sim_us('%t_y%',policy_bau))-1;
        yxf_pchg(r,s,policy_pchg,policy_bau)$(yxf_sim(r,s,'%t_y%',policy_bau) <> 0)       = (yxf_sim(r,s,'%t_y%',policy_pchg)/yxf_sim(r,s,'%t_y%',policy_bau))-1;
        yxf_pchg_r(r,policy_pchg,policy_bau)$(yxf_sim_r(r,'%t_y%',policy_bau) <> 0)       = (yxf_sim_r(r,'%t_y%',policy_pchg)/yxf_sim_r(r,'%t_y%',policy_bau))-1;
        yxf_pchg_us(policy_pchg,policy_bau)$(yxf_sim_us('%t_y%',policy_bau) <> 0)         = (yxf_sim_us('%t_y%',policy_pchg)/yxf_sim_us('%t_y%',policy_bau))-1;
        yt_pchg(r,s,policy_pchg,policy_bau)$(yt_sim(r,s,'%t_y%',policy_bau) <> 0)         = (yt_sim(r,s,'%t_y%',policy_pchg)/yt_sim(r,s,'%t_y%',policy_bau))-1;
        yt_pchg_r(r,policy_pchg,policy_bau)$(yt_sim_r(r,'%t_y%',policy_bau) <> 0)         = (yt_sim_r(r,'%t_y%',policy_pchg)/yt_sim_r(r,'%t_y%',policy_bau))-1;
        yt_pchg_us(policy_pchg,policy_bau)$(yt_sim_us('%t_y%',policy_bau) <> 0)           = (yt_sim_us('%t_y%',policy_pchg)/yt_sim_us('%t_y%',policy_bau))-1;
        a_pchg(r,s,policy_pchg,policy_bau)$(a_sim(r,s,'%t_y%',policy_bau) <> 0)           = (a_sim(r,s,'%t_y%',policy_pchg)/a_sim(r,s,'%t_y%',policy_bau))-1;

*       Activity by aggregate sector
        yd_pchg_sa(r,s_a,policy_pchg,policy_bau)$(yd_sim_sa(r,s_a,'%t_y%',policy_bau) <> 0)   = (yd_sim_sa(r,s_a,'%t_y%',policy_pchg)/yd_sim_sa(r,s_a,'%t_y%',policy_bau))-1;
        yxd_pchg_sa(r,s_a,policy_pchg,policy_bau)$(yxd_sim_sa(r,s_a,'%t_y%',policy_bau) <> 0) = (yxd_sim_sa(r,s_a,'%t_y%',policy_pchg)/yxd_sim_sa(r,s_a,'%t_y%',policy_bau))-1;
        yxf_pchg_sa(r,s_a,policy_pchg,policy_bau)$(yxf_sim_sa(r,s_a,'%t_y%',policy_bau) <> 0) = (yxf_sim_sa(r,s_a,'%t_y%',policy_pchg)/yxf_sim_sa(r,s_a,'%t_y%',policy_bau))-1;
        yt_pchg_sa(r,s_a,policy_pchg,policy_bau)$(yt_sim_sa(r,s_a,'%t_y%',policy_bau) <> 0)   = (yt_sim_sa(r,s_a,'%t_y%',policy_pchg)/yt_sim_sa(r,s_a,'%t_y%',policy_bau))-1;

*       Demand
        kd_pchg(r,s,policy_pchg,policy_bau)$(kd_sim(r,s,'%t_y%',policy_bau) <> 0)         = (kd_sim(r,s,'%t_y%',policy_pchg)/kd_sim(r,s,'%t_y%',policy_bau))-1;
        kd_pchg_r(r,policy_pchg,policy_bau)$(kd_sim_r(r,'%t_y%',policy_bau) <> 0)         = (kd_sim_r(r,'%t_y%',policy_pchg)/kd_sim_r(r,'%t_y%',policy_bau))-1;
        kd_pchg_us(policy_pchg,policy_bau)$(kd_sim_us('%t_y%',policy_bau) <> 0)           = (kd_sim_us('%t_y%',policy_pchg)/kd_sim_us('%t_y%',policy_bau))-1;
        ld_pchg(r,s,policy_pchg,policy_bau)$(ld_sim(r,s,'%t_y%',policy_bau) <> 0)         = (ld_sim(r,s,'%t_y%',policy_pchg)/ld_sim(r,s,'%t_y%',policy_bau))-1;
        ld_pchg_r(r,policy_pchg,policy_bau)$(ld_sim_r(r,'%t_y%',policy_bau) <> 0)         = (ld_sim_r(r,'%t_y%',policy_pchg)/ld_sim_r(r,'%t_y%',policy_bau))-1;
        ld_pchg_us(policy_pchg,policy_bau)$(ld_sim_us('%t_y%',policy_bau) <> 0)           = (ld_sim_us('%t_y%',policy_pchg)/ld_sim_us('%t_y%',policy_bau))-1;
        btd_pchg(r,s,policy_pchg,policy_bau)$(btd_sim(r,s,'%t_y%',policy_bau) <> 0)       = (btd_sim(r,s,'%t_y%',policy_pchg)/btd_sim(r,s,'%t_y%',policy_bau))-1;
        btd_pchg_r(r,policy_pchg,policy_bau)$(btd_sim_r(r,'%t_y%',policy_bau) <> 0)       = (btd_sim_r(r,'%t_y%',policy_pchg)/btd_sim_r(r,'%t_y%',policy_bau))-1;
        btd_pchg_us(policy_pchg,policy_bau)$(btd_sim_us('%t_y%',policy_bau) <> 0)         = (btd_sim_us('%t_y%',policy_pchg)/btd_sim_us('%t_y%',policy_bau))-1;
        id_pchg(r,g,s,policy_pchg,policy_bau)$(id_sim(r,g,s,'%t_y%',policy_bau) <> 0)     = (id_sim(r,g,s,'%t_y%',policy_pchg)/id_sim(r,g,s,'%t_y%',policy_bau))-1;
        id_pchg_g(g,policy_pchg,policy_bau)$(id_sim_g(g,'%t_y%',policy_bau) <> 0)         = (id_sim_g(g,'%t_y%',policy_pchg)/id_sim_g(g,'%t_y%',policy_bau))-1;
        id_pchg_r(r,policy_pchg,policy_bau)$(id_sim_r(r,'%t_y%',policy_bau) <> 0)         = (id_sim_r(r,'%t_y%',policy_pchg)/id_sim_r(r,'%t_y%',policy_bau))-1;
        id_pchg_us(policy_pchg,policy_bau)$(id_sim_us('%t_y%',policy_bau) <> 0)           = (id_sim_us('%t_y%',policy_pchg)/id_sim_us('%t_y%',policy_bau))-1;

        md_pchg(r,s,policy_pchg,policy_bau)$(md_sim(r,s,'%t_y%',policy_bau) <> 0)         = (md_sim(r,s,'%t_y%',policy_pchg)/md_sim(r,s,'%t_y%',policy_bau))-1;
        md_pchg_r(r,policy_pchg,policy_bau)$(md_sim_r(r,'%t_y%',policy_bau) <> 0)         = (md_sim_r(r,'%t_y%',policy_pchg)/md_sim_r(r,'%t_y%',policy_bau))-1;
        md_pchg_us(policy_pchg,policy_bau)$(md_sim_us('%t_y%',policy_bau) <> 0)           = (md_sim_us('%t_y%',policy_pchg)/md_sim_us('%t_y%',policy_bau))-1;
        mf_pchg(r,s,policy_pchg,policy_bau)$(mf_sim(r,s,'%t_y%',policy_bau) <> 0)         = (mf_sim(r,s,'%t_y%',policy_pchg)/mf_sim(r,s,'%t_y%',policy_bau))-1;
        mf_pchg_r(r,policy_pchg,policy_bau)$(mf_sim_r(r,'%t_y%',policy_bau) <> 0)         = (mf_sim_r(r,'%t_y%',policy_pchg)/mf_sim_r(r,'%t_y%',policy_bau))-1;
        mf_pchg_us(policy_pchg,policy_bau)$(mf_sim_us('%t_y%',policy_bau) <> 0)           = (mf_sim_us('%t_y%',policy_pchg)/mf_sim_us('%t_y%',policy_bau))-1;

        cd_pchg(r,h,policy_pchg,policy_bau)$(cd_sim(r,h,'%t_y%',policy_bau) <> 0)         = (cd_sim(r,h,'%t_y%',policy_pchg)/cd_sim(r,h,'%t_y%',policy_bau))-1;
        cd_pchg_h(h,policy_pchg,policy_bau)$(cd_sim_h(h,'%t_y%',policy_bau) <> 0)         = (cd_sim_h(h,'%t_y%',policy_pchg)/cd_sim_h(h,'%t_y%',policy_bau))-1;
        cd_pchg_r(r,policy_pchg,policy_bau)$(cd_sim_r(r,'%t_y%',policy_bau) <> 0)         = (cd_sim_r(r,'%t_y%',policy_pchg)/cd_sim_r(r,'%t_y%',policy_bau))-1;
        cd_pchg_us(policy_pchg,policy_bau)$(cd_sim_us('%t_y%',policy_bau) <> 0)           = (cd_sim_us('%t_y%',policy_pchg)/cd_sim_us('%t_y%',policy_bau))-1;
        i_pchg_r(r,policy_pchg,policy_bau)$(i_sim_r(r,'%t_y%',policy_bau) <> 0)           = (i_sim_r(r,'%t_y%',policy_pchg)/i_sim_r(r,'%t_y%',policy_bau))-1;
        i_pchg_us(policy_pchg,policy_bau)$(i_sim_us('%t_y%',policy_bau) <> 0)             = (i_sim_us('%t_y%',policy_pchg)/i_sim_us('%t_y%',policy_bau))-1;
        g_pchg(r,pub,policy_pchg,policy_bau)$(g_sim(r,pub,'%t_y%',policy_bau) <> 0)       = (g_sim(r,pub,'%t_y%',policy_pchg)/g_sim(r,pub,'%t_y%',policy_bau))-1;
        g_pchg_r(r,policy_pchg,policy_bau)$(g_sim_r(r,'%t_y%',policy_bau) <> 0)           = (g_sim_r(r,'%t_y%',policy_pchg)/g_sim_r(r,'%t_y%',policy_bau))-1;
        g_pchg_us(policy_pchg,policy_bau)$(g_sim_us('%t_y%',policy_bau) <> 0)             = (g_sim_us('%t_y%',policy_pchg)/g_sim_us('%t_y%',policy_bau))-1;

*       Additional
        trev_pchg_r(r,policy_pchg,policy_bau)$(trev_sim_r(r,'%t_y%',policy_bau) <> 0)     = (trev_sim_r(r,'%t_y%',policy_pchg)/trev_sim_r(r,'%t_y%',policy_bau))-1;
        trev_pchg_us(policy_pchg,policy_bau)$(trev_sim_us('%t_y%',policy_bau) <> 0)       = (trev_sim_us('%t_y%',policy_pchg)/trev_sim_us('%t_y%',policy_bau))-1;
        grpfc_pchg_r(r,policy_pchg,policy_bau)$(grpfc_sim_r(r,'%t_y%',policy_bau) <> 0)   = (grpfc_sim_r(r,'%t_y%',policy_pchg)/grpfc_sim_r(r,'%t_y%',policy_bau))-1;
        grpfc_pchg_us(policy_pchg,policy_bau)$(grpfc_sim_us('%t_y%',policy_bau) <> 0)     = (grpfc_sim_us('%t_y%',policy_pchg)/grpfc_sim_us('%t_y%',policy_bau))-1;
        grpva_pchg_r(r,policy_pchg,policy_bau)$(grpva_sim_r(r,'%t_y%',policy_bau) <> 0)   = (grpva_sim_r(r,'%t_y%',policy_pchg)/grpva_sim_r(r,'%t_y%',policy_bau))-1;
        grpva_pchg_us(policy_pchg,policy_bau)$(grpva_sim_us('%t_y%',policy_bau) <> 0)     = (grpva_sim_us('%t_y%',policy_pchg)/grpva_sim_us('%t_y%',policy_bau))-1;
    );
);


* ---------------------------------------------------------------------------- *
*        Display and Save Results:
* ---------------------------------------------------------------------------- *

* Model parameters

$label display_recursive

$if %t_y%==2010 $goto display_psm

DISPLAY
        simtp, ts, tts,
        kstock_sim, lstock_sim, evok_sim, evol_sim, incadj_cal_sim,
        incadj_sim, rhscale_sim, ca_sim, hinv_sim, evolp_sim, evokp_sim,
        lpidx_sim
;

* PSM parameters

$label display_psm

$if not %s_pw%==1 $goto display_drm

DISPLAY
      psmtfp_sim, psmtfp_fa_sim
;

* Demand response parameters
$label display_drm
$if not %s_dr%==1 $goto display_results

DISPLAY
      shk_dr_pchg_prod_sim, shk_dr_pchg_consum_sim,
      dre_prod_sim, dre_hh_sim
;

*        Simulation level output

$label display_results

DISPLAY
        drate, irate, kprate, lprate, poprate,
        py_sim, pa_sim, pc_sim, pn_sim, pinv_sim,
        pgov_sim, pk_sim, pl_sim, pbt_sim,
        yd_sim, yd_sim_r, yd_sim_us,
        vyd_sim,
        yxd_sim, yxd_sim_r, yxd_sim_us,
        yxf_sim, yxf_sim_r, yxf_sim_us,
        yt_sim, yt_sim_r, yt_sim_us,
        yd_sim_sa, yxd_sim_sa, yxf_sim_sa, yt_sim_sa,
        yt_share_sa, yt_share_r,
        kd_sim, kd_sim_r, kd_sim_us,
        ld_sim, ld_sim_r, ld_sim_us,
        btd_sim, btd_sim_r, btd_sim_us,
        id_sim, id_sim_g, id_sim_r, id_sim_us,
        trev_sim_r,
        cd_sim, cd_sim_r, cd_sim_us,
        a_sim,
        md_sim, md_sim_r, md_sim_us,
        mf_sim, mf_sim_r, mf_sim_us,
        cd_sim, cd_sim_h, cd_sim_r, cd_sim_us,
        g_sim, g_sim_r, g_sim_us,
        i_sim_r, i_sim_us,
        trev_sim_r, trev_sim_us,
        grpfc_sim_r, grpfc_sim_us,
        grpva_sim_r, grpva_sim_us
;

execute_unload '%save_valuelevel%'
        forloop, policy, policy_bau,
        zones, econr, zones_map, econr_map,
        r_psm_zones, r_psm_econr,
        lpidx_sim, kstock_sim, lstock_sim,
        py_sim, pa_sim, pc_sim, pn_sim, pinv_sim,
        pgov_sim, pk_sim, pl_sim, pbt_sim,
        yd_sim, yd_sim_r, yd_sim_us,
        vyd_sim,
        yxd_sim, yxd_sim_r, yxd_sim_us,
        yxf_sim, yxf_sim_r, yxf_sim_us,
        yt_sim, yt_sim_r, yt_sim_us,
        yt_share_sa, yt_share_r,
        yt_share_sa, yt_share_r,
        yd_sim_sa, yxd_sim_sa, yxf_sim_sa, yt_sim_sa,
        a_sim,
        kd_sim, kd_sim_r, kd_sim_us,
        ld_sim, ld_sim_r, ld_sim_us,
        btd_sim, btd_sim_r, btd_sim_us,
        id_sim, id_sim_g, id_sim_r, id_sim_us,
        cd_sim, cd_sim_r, cd_sim_us,
        md_sim, md_sim_r, md_sim_us,
        mf_sim, mf_sim_r, mf_sim_us,
        cd_sim, cd_sim_h, cd_sim_r, cd_sim_us,
        g_sim, g_sim_r, g_sim_us,
        i_sim_r, i_sim_us,
        trev_sim_r, trev_sim_us,
        grpfc_sim_r, grpfc_sim_us,
        grpva_sim_r, grpva_sim_us
;

*        Simulation percentage change output
DISPLAY
        py_pchg, pa_pchg, pc_pchg, pn_pchg, pinv_pchg,
        pgov_pchg, pk_pchg, pl_pchg, pbt_pchg,
        yd_pchg, yd_pchg_r, yd_pchg_us,
        vyd_pchg,
        yxd_pchg, yxd_pchg_r, yxd_pchg_us,
        yxf_pchg, yxf_pchg_r, yxf_pchg_us,
        yt_pchg, yt_pchg_r, yt_pchg_us,
        yd_pchg_sa, yxd_pchg_sa, yxf_pchg_sa, yt_pchg_sa,
        a_pchg,
        kd_pchg, kd_pchg_r, kd_pchg_us,
        ld_pchg, ld_pchg_r, ld_pchg_us,
        btd_pchg, btd_pchg_r, btd_pchg_us,
        id_pchg, id_pchg_g, id_pchg_r, id_pchg_us,
        md_pchg, md_pchg_r, md_pchg_us,
        mf_pchg, mf_pchg_r, mf_pchg_us,
        cd_pchg, cd_pchg_h, cd_pchg_r, cd_pchg_us,
        i_pchg_r, i_pchg_us,
        g_pchg, g_pchg_r, g_pchg_us,
        trev_pchg_r, trev_pchg_us,
        grpfc_pchg_r, grpfc_pchg_us,
        grpva_pchg_r, grpva_pchg_us
;

execute_unload '%save_valuepchg%'
        forloop, policy, policy_bau, policy_pchg,
        zones, econr, zones_map, econr_map,
        r_psm_zones, r_psm_econr,
        py_pchg, pa_pchg, pc_pchg, pn_pchg, pinv_pchg,
        pgov_pchg, pk_pchg, pl_pchg, pbt_pchg,
        yd_pchg, yd_pchg_r, yd_pchg_us,
        vyd_pchg,
        yxd_pchg, yxd_pchg_r, yxd_pchg_us,
        yxf_pchg, yxf_pchg_r, yxf_pchg_us,
        yt_pchg, yt_pchg_r, yt_pchg_us,
        yd_pchg_sa, yxd_pchg_sa, yxf_pchg_sa, yt_pchg_sa,
        a_pchg,
        kd_pchg, kd_pchg_r, kd_pchg_us,
        ld_pchg, ld_pchg_r, ld_pchg_us,
        btd_pchg, btd_pchg_r, btd_pchg_us,
        id_pchg, id_pchg_g, id_pchg_r, id_pchg_us,
        md_pchg, md_pchg_r, md_pchg_us,
        mf_pchg, mf_pchg_r, mf_pchg_us,
        cd_pchg, cd_pchg_h, cd_pchg_r, cd_pchg_us,
        i_pchg_r, i_pchg_us,
        g_pchg, g_pchg_r, g_pchg_us,
        trev_pchg_r, trev_pchg_us,
        grpfc_pchg_r, grpfc_pchg_us,
        grpva_pchg_r, grpva_pchg_us
;

*        Simulation growth rate output

$label display_cagr

$if %t_y%==2010 $goto save_model_parameters

DISPLAY
        yd_cagr,  yd_cagr_r,  yd_cagr_us,
        yxd_cagr, yxd_cagr_r, yxd_cagr_us,
        yxf_cagr, yxf_cagr_r, yxf_cagr_us,
        yt_cagr,  yt_cagr_r,  yt_cagr_us,
        yd_cagr_sa, yxd_cagr_sa, yxf_cagr_sa, yt_cagr_sa,
        grpfc_cagr_r, grpfc_cagr_us,
        grpva_cagr_r, grpva_cagr_us,
        yd_tcagr, yd_tcagr_r,  yd_tcagr_us,
        yxd_tcagr, yxd_tcagr_r, yxd_tcagr_us,
        yxf_tcagr, yxf_tcagr_r, yxf_tcagr_us,
        yd_tcagr_sa, yxd_tcagr_sa, yxf_tcagr_sa, yt_tcagr_sa,
        grpfc_tcagr_r, grpfc_tcagr_us,
        grpva_tcagr_r, grpva_tcagr_us
;

execute_unload '%save_growthrate%'
        forloop, policy, policy_bau,
        zones, econr, zones_map, econr_map,
        r_psm_zones, r_psm_econr,
        yd_cagr, yd_cagr_r, yd_cagr_us,
        yxd_cagr, yxd_cagr_r, yxd_cagr_us,
        yxf_cagr, yxf_cagr_r, yxf_cagr_us,
        yt_cagr, yt_cagr_r, yt_cagr_us,
        yd_cagr_sa, yxd_cagr_sa, yxf_cagr_sa, yt_cagr_sa,
        grpfc_cagr_r, grpfc_cagr_us,
        grpva_cagr_r, grpva_cagr_us,
        yd_tcagr, yd_tcagr_r, yd_tcagr_us,
        yxd_tcagr, yxd_tcagr_r, yxd_tcagr_us,
        yxf_tcagr, yxf_tcagr_r, yxf_tcagr_us,
        yt_tcagr, yt_tcagr_r, yt_tcagr_us,
        yd_tcagr_sa, yxd_tcagr_sa, yxf_tcagr_sa, yt_tcagr_sa,
        grpfc_tcagr_r, grpfc_tcagr_us,
        grpva_tcagr_r, grpva_tcagr_us
;

* ---------------------------------------------------------------------------- *
*        Save File with model parameter assumptions:
* ---------------------------------------------------------------------------- *

$label save_model_parameters

file results /'%save_parameters%'/;
results.pc = 5;
results.nd = 6;
put results;
* header
put "parameter", "title", "notes" "value"/;
* model parameters
put "target", "IMPLAN Aggregation Dataset",        "implan440-30sector-12region", "%target%"/;
put "t_r",    "Number of Regions",                 "default == 12",                  %t_r%/;
put "t_s",    "Number of Sectors",                 "default == 30",                  %t_s%/;
put "t_y",    "The Final Time Period",             "default == 2050",                %t_y%/;
put "t_p",    "Shared Social Pathway Model",       "default == SSP2",           "SSP%t_p%"/;
put "s_wd",   "WDRG coupling",                     "== 1 on; == 0 off",             %s_wd%/;
put "s_pw",   "PSM w/o Sullivan coupling",         "== 1 on; == 0 off",             %s_pw%/;
put "s_pf",   "PSM w/Sullivan coupling",           "== 1 on; == 0 off",             %s_pf%/;
put "s_dr",   "Demand Response Emulator coupling", "== 1 on; == 0 off",             %s_dr%/;
* Macro parameters
put "drate",  "Depreciation rate",                 "default == 0.03",               drate/;
put "irate",  "Interest rate",                     "default == 0.05",               irate/;
put "kprate", "National capital productivity",     "default == 0.11",               kprate/;
put "lprate", "National labor productivity",       "default == 0.016",              lprate/;

putclose;

*display py_sim.tl;

file results1 /'%save1_parameters%'/;
results1.pc = 5;
results1.nd = 6;
put results1;
* header
put                               "parameter",   "parameter.title",  "simulation", "year", "region", "sector.aggregation", "sector", "household", "government", "from", "value" /;
* capital and labor stock
*loop((forloop,simtp,r),     put   "kstock_sim",    kstock_sim.ts,    forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",    "",      "",    kstock_sim(r,simtp,forloop)    / );
*loop((forloop,simtp,r),     put   "lstock_sim",    lstock_sim.ts,    forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",    "",      "",    lstock_sim(r,simtp,forloop)    / );
* relative prices
loop((forloop,simtp,r,s),   put   "py_sim",        py_sim.ts,        forloop.tl,  simtp.tl,  r.tl,  "",      s.tl,  "",    "",      "",    py_sim(r,s,simtp,forloop)      / );
loop((forloop,simtp,r,s),   put   "pa_sim",        pa_sim.ts,        forloop.tl,  simtp.tl,  r.tl,  "",      s.tl,  "",    "",      "",    pa_sim(r,s,simtp,forloop)      / );
loop((forloop,simtp,r,h),   put   "pc_sim",        pc_sim.ts,        forloop.tl,  simtp.tl,  r.tl,  "",      "",    h.tl,  "",      "",    pc_sim(r,h,simtp,forloop)      / );
loop((forloop,simtp,s),     put   "pn_sim",        pn_sim.ts,        forloop.tl,  simtp.tl,  "",    "",      s.tl,  "",    "",      "",    pn_sim(s,simtp,forloop)        / );
loop((forloop,simtp,r),     put   "pinv_sim",      pinv_sim.ts,      forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",    "",      "",    pinv_sim(r,simtp,forloop)      / );
loop((forloop,simtp,r,pub), put   "pgov_sim",      pgov_sim.ts,      forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",    pub.tl,  "",    pgov_sim(r,pub,simtp,forloop)  / );
loop((forloop,simtp,r),     put   "pk_sim",        pk_sim.ts,        forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",    "",      "",    pk_sim(r,simtp,forloop)        / );
loop((forloop,simtp,r),     put   "pl_sim",        pl_sim.ts,        forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",    "",      "",    pl_sim(r,simtp,forloop)        / );
loop((forloop,simtp,r),     put   "pbt_sim",       pbt_sim.ts,       forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",    "",      "",    pbt_sim(r,simtp,forloop)       / );
* activity - production
loop((forloop,simtp,r,s),   put   "yd_sim",        yd_sim.ts,        forloop.tl,  simtp.tl,  r.tl,  "",      s.tl,  "",    "",      "",    yd_sim(r,s,simtp,forloop)      / );
loop((forloop,simtp,r),     put   "yd_sim_r",      yd_sim_r.ts,      forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",    "",      "",    yd_sim_r(r,simtp,forloop)      / );
loop((forloop,simtp),       put   "yd_sim_us",     yd_sim_us.ts,     forloop.tl,  simtp.tl,  "",    "",      "",    "",    "",      "",    yd_sim_us(simtp,forloop)       / );
loop((forloop,simtp,r,s),   put   "vyd_sim",       vyd_sim.ts,       forloop.tl,  simtp.tl,  r.tl,  "",      s.tl,  "",    "",      "",    vyd_sim(r,s,simtp,forloop)     / );
loop((forloop,simtp,r,s),   put   "yxd_sim",       yxd_sim.ts,       forloop.tl,  simtp.tl,  r.tl,  "",      s.tl,  "",    "",      "",    yxd_sim(r,s,simtp,forloop)     / );
loop((forloop,simtp,r),     put   "yxd_sim_r",     yxd_sim_r.ts,     forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",    "",      "",    yxd_sim_r(r,simtp,forloop)     / );
loop((forloop,simtp),       put   "yxd_sim_us",    yxd_sim_us.ts,    forloop.tl,  simtp.tl,  "",    "",      "",    "",    "",      "",    yxd_sim_us(simtp,forloop)      / );
loop((forloop,simtp,r,s),   put   "yxf_sim",       yxf_sim.ts,       forloop.tl,  simtp.tl,  r.tl,  "",      s.tl,  "",    "",      "",    yxf_sim(r,s,simtp,forloop)     / );
loop((forloop,simtp,r),     put   "yxf_sim_r",     yxf_sim_r.ts,     forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",    "",      "",    yxf_sim_r(r,simtp,forloop)     / );
loop((forloop,simtp),       put   "yxf_sim_us",    yxf_sim_us.ts,    forloop.tl,  simtp.tl,  "",    "",      "",    "",    "",      "",    yxf_sim_us(simtp,forloop)      / );
loop((forloop,simtp,r,s),   put   "yt_sim",        yt_sim.ts,        forloop.tl,  simtp.tl,  r.tl,  "",      s.tl,  "",    "",      "",    yt_sim(r,s,simtp,forloop)      / );
loop((forloop,simtp,r),     put   "yt_sim_r",      yt_sim_r.ts,      forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",    "",      "",    yt_sim_r(r,simtp,forloop)      / );
loop((forloop,simtp),       put   "yt_sim_us",     yt_sim_us.ts,     forloop.tl,  simtp.tl,  "",    "",      "",    "",    "",      "",    yt_sim_us(simtp,forloop)       / );
loop((forloop,simtp,r,s),   put   "a_sim",         a_sim.ts,         forloop.tl,  simtp.tl,  r.tl,  "",      s.tl,  "",    "",      "",    a_sim(r,s,simtp,forloop)       / );
* aggregate production
loop((forloop,simtp,r,s_a), put   "yd_sim_sa",     yd_sim_sa.ts,     forloop.tl,  simtp.tl,  r.tl,  s_a.tl,  "",    "",    "",      "",    yd_sim_sa(r,s_a,simtp,forloop) / );
loop((forloop,simtp,r,s_a), put   "yxd_sim_sa",    yxd_sim_sa.ts,    forloop.tl,  simtp.tl,  r.tl,  s_a.tl,  "",    "",    "",      "",    yxd_sim_sa(r,s_a,simtp,forloop) / );
loop((forloop,simtp,r,s_a), put   "yxf_sim_sa",    yxf_sim_sa.ts,    forloop.tl,  simtp.tl,  r.tl,  s_a.tl,  "",    "",    "",      "",    yxf_sim_sa(r,s_a,simtp,forloop) / );
loop((forloop,simtp,r,s_a), put   "yt_sim_sa",     yt_sim_sa.ts,     forloop.tl,  simtp.tl,  r.tl,  s_a.tl,  "",    "",    "",      "",    yt_sim_sa(r,s_a,simtp,forloop) / );
* demand - factor
loop((forloop,simtp,r,s),   put   "kd_sim",        kd_sim.ts,        forloop.tl,  simtp.tl,  r.tl, "",       s.tl,  "",    "",      "",    kd_sim(r,s,simtp,forloop) / );
loop((forloop,simtp,r),     put   "kd_sim_r",      kd_sim_r.ts,      forloop.tl,  simtp.tl,  r.tl, "",       "",    "",    "",      "",    kd_sim_r(r,simtp,forloop) / );
loop((forloop,simtp),       put   "kd_sim_us",     kd_sim_us.ts,     forloop.tl,  simtp.tl,  "",   "",       "",    "",    "",      "",    kd_sim_us(simtp,forloop) / );
loop((forloop,simtp,r,s),   put   "ld_sim",        ld_sim.ts,        forloop.tl,  simtp.tl,  r.tl, "",       s.tl,  "",    "",      "",    ld_sim(r,s,simtp,forloop) / );
loop((forloop,simtp,r),     put   "ld_sim_r",      ld_sim_r.ts,      forloop.tl,  simtp.tl,  r.tl, "",       "",    "",    "",      "",    ld_sim_r(r,simtp,forloop) / );
loop((forloop,simtp),       put   "ld_sim_us",     ld_sim_us.ts,     forloop.tl,  simtp.tl,  "",   "",       "",    "",    "",      "",    ld_sim_us(simtp,forloop) / );
loop((forloop,simtp,r,s),   put   "btd_sim",       btd_sim.ts,       forloop.tl,  simtp.tl,  r.tl, "",       s.tl,  "",    "",      "",    btd_sim(r,s,simtp,forloop) / );
loop((forloop,simtp,r),     put   "btd_sim_r",     btd_sim_r.ts,     forloop.tl,  simtp.tl,  r.tl, "",       "",    "",    "",      "",    btd_sim_r(r,simtp,forloop) / );
loop((forloop,simtp),       put   "btd_sim_us",    btd_sim_us.ts,    forloop.tl,  simtp.tl,  "",   "",       "",    "",    "",      "",    btd_sim_us(simtp,forloop) / );
* demand - intermediate
loop((forloop,simtp,r,s,g), put   "id_sim",        id_sim.ts,        forloop.tl,  simtp.tl,  r.tl, "",       s.tl,  "",    "",      g.tl,  id_sim(r,g,s,simtp,forloop) / );
loop((forloop,simtp,g),     put   "id_sim_g",      id_sim_g.ts,      forloop.tl,  simtp.tl,  "",   "",       "",    "",    "",      g.tl,  id_sim_g(g,simtp,forloop) / );
loop((forloop,simtp,r),     put   "id_sim_r",      id_sim_r.ts,      forloop.tl,  simtp.tl,  r.tl, "",       "",    "",    "",      "",    id_sim_r(r,simtp,forloop) / );
loop((forloop,simtp),       put   "id_sim_us",     id_sim_us.ts,     forloop.tl,  simtp.tl,  "",   "",       "",    "",    "",      "",    id_sim_us(simtp,forloop) / );
* demand - export and import
loop((forloop,simtp,r,s),   put   "md_sim",        md_sim.ts,        forloop.tl,  simtp.tl,  r.tl, "",       s.tl,  "",    "",      "",    md_sim(r,s,simtp,forloop) / );
loop((forloop,simtp,r),     put   "md_sim_r",      md_sim_r.ts,      forloop.tl,  simtp.tl,  r.tl, "",       "",    "",    "",      "",    md_sim_r(r,simtp,forloop) / );
loop((forloop,simtp),       put   "md_sim_us",     md_sim_us.ts,     forloop.tl,  simtp.tl,  "",   "",       "",    "",    "",      "",    md_sim_us(simtp,forloop) / );
loop((forloop,simtp,r,s),   put   "mf_sim",        mf_sim.ts,        forloop.tl,  simtp.tl,  r.tl, "",       s.tl,  "",    "",      "",    mf_sim(r,s,simtp,forloop) / );
loop((forloop,simtp,r),     put   "mf_sim_r",      mf_sim_r.ts,      forloop.tl,  simtp.tl,  r.tl, "",       "",    "",    "",      "",    mf_sim_r(r,simtp,forloop) / );
loop((forloop,simtp),       put   "mf_sim_us",     mf_sim_us.ts,     forloop.tl,  simtp.tl,  "",   "",       "",    "",    "",      "",    mf_sim_us(simtp,forloop) / );
* demand - consumption, investiment, and government
loop((forloop,simtp,r,h),   put   "cd_sim",        cd_sim.ts,        forloop.tl,  simtp.tl,  r.tl, "",       "",    h.tl,  "",      "",    cd_sim(r,h,simtp,forloop) / );
loop((forloop,simtp,h),     put   "cd_sim_h",      cd_sim_h.ts,      forloop.tl,  simtp.tl,  "",   "",       "",    h.tl,  "",      "",    cd_sim_h(h,simtp,forloop) / );
loop((forloop,simtp,r),     put   "cd_sim_r",      cd_sim_r.ts,      forloop.tl,  simtp.tl,  r.tl, "",       "",    "",    "",      "",    cd_sim_r(r,simtp,forloop) / );
loop((forloop,simtp),       put   "cd_sim_us",     cd_sim_us.ts,     forloop.tl,  simtp.tl,  "",   "",       "",    "",    "",      "",    cd_sim_us(simtp,forloop) / );
loop((forloop,simtp,r),     put   "i_sim_r",       i_sim_r.ts,       forloop.tl,  simtp.tl,  r.tl, "",       "",    "",    "",      "",    i_sim_r(r,simtp,forloop) / );
loop((forloop,simtp),       put   "i_sim_us",      i_sim_us.ts,      forloop.tl,  simtp.tl,  "",   "",       "",    "",    "",      "",    i_sim_us(simtp,forloop) / );
loop((forloop,simtp,r,pub), put   "g_sim",         g_sim.ts,         forloop.tl,  simtp.tl,  r.tl, "",       "",    "",    pub.tl,  "",    g_sim(r,pub,simtp,forloop) / );
loop((forloop,simtp,r),     put   "g_sim_r",       g_sim_r.ts,       forloop.tl,  simtp.tl,  r.tl, "",       "",    "",    "",      "",    g_sim_r(r,simtp,forloop) / );
loop((forloop,simtp),       put   "g_sim_us",      g_sim_us.ts,      forloop.tl,  simtp.tl,  "",   "",       "",    "",    "",      "",    g_sim_us(simtp,forloop) / );
* tax revenue
loop((forloop,simtp,r),     put   "trev_sim_r",    trev_sim_r.ts,    forloop.tl,  simtp.tl,  r.tl, "",       "",    "",    "",      "",    trev_sim_r(r,simtp,forloop) / );
loop((forloop,simtp),       put   "trev_sim_us",   trev_sim_us.ts,   forloop.tl,  simtp.tl,  "",   "",       "",    "",    "",      "",    trev_sim_us(simtp,forloop) / );
* gross regional product
loop((forloop,simtp,r),     put   "grpfc_sim_r",   grpfc_sim_r.ts,   forloop.tl,  simtp.tl, r.tl,  "",       "",    "",    "",      "",    grpfc_sim_r(r,simtp,forloop) / );
loop((forloop,simtp),       put   "grpfc_sim_us",  grpfc_sim_us.ts,  forloop.tl,  simtp.tl, "",    "",       "",    "",    "",      "",    grpfc_sim_us(simtp,forloop) / );
loop((forloop,simtp,r),     put   "grpva_sim_r",   grpva_sim_r.ts,   forloop.tl,  simtp.tl, r.tl,  "",       "",    "",    "",      "",    grpva_sim_r(r,simtp,forloop) / );
loop((forloop,simtp),       put   "grpva_sim_us",  grpva_sim_us.ts,  forloop.tl,  simtp.tl, "",    "",       "",    "",    "",      "",    grpva_sim_us(simtp,forloop) / );

putclose;

file results2 /'%save2_parameters%'/;
results2.pc = 5;
results2.nd = 6;
put results2;
* header
put                     "parameter",     "parameter.title",  "simulation",  "bau simulation", "region", "sector.aggregation", "sector", "household", "government", "from", "value" /;
* relative prices
loop((policy,policy_bau,r,s),    put  "py_pchg",        py_pchg.ts,        policy.tl,  policy_bau.tl,  r.tl,  "",      s.tl,  "",    "",      "",    py_pchg(r,s,policy,policy_bau)        / );
loop((policy,policy_bau,r,s),    put  "pa_pchg",        pa_pchg.ts,        policy.tl,  policy_bau.tl,  r.tl,  "",      s.tl,  "",    "",      "",    pa_pchg(r,s,policy,policy_bau)        / );
loop((policy,policy_bau,r,h),    put  "pc_pchg",        pc_pchg.ts,        policy.tl,  policy_bau.tl,  r.tl,  "",      "",    h.tl,  "",      "",    pc_pchg(r,h,policy,policy_bau)        / );
loop((policy,policy_bau,s),      put  "pn_pchg",        pn_pchg.ts,        policy.tl,  policy_bau.tl,  "",    "",      s.tl,  "",    "",      "",    pn_pchg(s,policy,policy_bau)          / );
loop((policy,policy_bau,r),      put  "pinv_pchg",      pinv_pchg.ts,      policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    pinv_pchg(r,policy,policy_bau)        / );
loop((policy,policy_bau,r,pub),  put  "pgov_pchg",      pgov_pchg.ts,      policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    pub.tl,  "",    pgov_pchg(r,pub,policy,policy_bau)    / );
loop((policy,policy_bau,r),      put  "pk_pchg",        pk_pchg.ts,        policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    pk_pchg(r,policy,policy_bau)          / );
loop((policy,policy_bau,r),      put  "pl_pchg",        pl_pchg.ts,        policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    pl_pchg(r,policy,policy_bau)          / );
loop((policy,policy_bau,r),      put  "pbt_pchg",       pbt_pchg.ts,       policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    pbt_pchg(r,policy,policy_bau)         / );
* activity - production
loop((policy,policy_bau,r,s),    put  "yd_pchg",        yd_pchg.ts,        policy.tl,  policy_bau.tl,  r.tl,  "",      s.tl,  "",    "",      "",    yd_pchg(r,s,policy,policy_bau)        / );
loop((policy,policy_bau,r),      put  "yd_pchg_r",      yd_pchg_r.ts,      policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    yd_pchg_r(r,policy,policy_bau)        / );
loop((policy,policy_bau),        put  "yd_pchg_us",     yd_pchg_us.ts,     policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    yd_pchg_us(policy,policy_bau)         / );
loop((policy,policy_bau,r,s),    put  "vyd_pchg",       vyd_pchg.ts,       policy.tl,  policy_bau.tl,  r.tl,  "",      s.tl,  "",    "",      "",    vyd_pchg(r,s,policy,policy_bau)       / );
loop((policy,policy_bau,r,s),    put  "yxd_pchg",       yxd_pchg.ts,       policy.tl,  policy_bau.tl,  r.tl,  "",      s.tl,  "",    "",      "",    yxd_pchg(r,s,policy,policy_bau)       / );
loop((policy,policy_bau,r),      put  "yxd_pchg_r",     yxd_pchg_r.ts,     policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    yxd_pchg_r(r,policy,policy_bau)       / );
loop((policy,policy_bau),        put  "yxd_pchg_us",    yxd_pchg_us.ts,    policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    yxd_pchg_us(policy,policy_bau)        / );
loop((policy,policy_bau,r,s),    put  "yxf_pchg",       yxf_pchg.ts,       policy.tl,  policy_bau.tl,  r.tl,  "",      s.tl,  "",    "",      "",    yxf_pchg(r,s,policy,policy_bau)       / );
loop((policy,policy_bau,r),      put  "yxf_pchg_r",     yxf_pchg_r.ts,     policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    yxf_pchg_r(r,policy,policy_bau)       / );
loop((policy,policy_bau),        put  "yxf_pchg_us",    yxf_pchg_us.ts,    policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    yxf_pchg_us(policy,policy_bau)        / );
loop((policy,policy_bau,r,s),    put  "yt_pchg",        yt_pchg.ts,        policy.tl,  policy_bau.tl,  r.tl,  "",      s.tl,  "",    "",      "",    yt_pchg(r,s,policy,policy_bau)        / );
loop((policy,policy_bau,r),      put  "yt_pchg_r",      yt_pchg_r.ts,      policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    yt_pchg_r(r,policy,policy_bau)        / );
loop((policy,policy_bau),        put  "yt_pchg_us",     yt_pchg_us.ts,     policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    yt_pchg_us(policy,policy_bau)         / );
loop((policy,policy_bau,r,s),    put  "a_pchg",         a_pchg.ts,         policy.tl,  policy_bau.tl,  r.tl,  "",      s.tl,  "",    "",      "",    a_pchg(r,s,policy,policy_bau)         / );
* aggregate production
loop((policy,policy_bau,r,s_a),  put  "yd_pchg_sa",     yd_pchg_sa.ts,     policy.tl,  policy_bau.tl,  r.tl,  s_a.tl,  "",    "",    "",      "",    yd_pchg_sa(r,s_a,policy,policy_bau)   / );
loop((policy,policy_bau,r,s_a),  put  "yxd_pchg_sa",    yxd_pchg_sa.ts,    policy.tl,  policy_bau.tl,  r.tl,  s_a.tl,  "",    "",    "",      "",    yxd_pchg_sa(r,s_a,policy,policy_bau)  / );
loop((policy,policy_bau,r,s_a),  put  "yxf_pchg_sa",    yxf_pchg_sa.ts,    policy.tl,  policy_bau.tl,  r.tl,  s_a.tl,  "",    "",    "",      "",    yxf_pchg_sa(r,s_a,policy,policy_bau)  / );
loop((policy,policy_bau,r,s_a),  put  "yt_pchg_sa",     yt_pchg_sa.ts,     policy.tl,  policy_bau.tl,  r.tl,  s_a.tl,  "",    "",    "",      "",    yt_pchg_sa(r,s_a,policy,policy_bau)   / );
* demand - factor
loop((policy,policy_bau,r,s),    put  "kd_pchg",        kd_pchg.ts,        policy.tl,  policy_bau.tl,  r.tl,  "",      s.tl,  "",    "",      "",    kd_pchg(r,s,policy,policy_bau)        / );
loop((policy,policy_bau,r),      put  "kd_pchg_r",      kd_pchg_r.ts,      policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    kd_pchg_r(r,policy,policy_bau)        / );
loop((policy,policy_bau),        put  "kd_pchg_us",     kd_pchg_us.ts,     policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    kd_pchg_us(policy,policy_bau)         / );
loop((policy,policy_bau,r,s),    put  "ld_pchg",        ld_pchg.ts,        policy.tl,  policy_bau.tl,  r.tl,  "",      s.tl,  "",    "",      "",    ld_pchg(r,s,policy,policy_bau)        / );
loop((policy,policy_bau,r),      put  "ld_pchg_r",      ld_pchg_r.ts,      policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    ld_pchg_r(r,policy,policy_bau)        / );
loop((policy,policy_bau),        put  "ld_pchg_us",     ld_pchg_us.ts,     policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    ld_pchg_us(policy,policy_bau)         / );
loop((policy,policy_bau,r,s),    put  "btd_pchg",       btd_pchg.ts,       policy.tl,  policy_bau.tl,  r.tl,  "",      s.tl,  "",    "",      "",    btd_pchg(r,s,policy,policy_bau)       / );
loop((policy,policy_bau,r),      put  "btd_pchg_r",     btd_pchg_r.ts,     policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    btd_pchg_r(r,policy,policy_bau)       / );
loop((policy,policy_bau),        put  "btd_pchg_us",    btd_pchg_us.ts,    policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    btd_pchg_us(policy,policy_bau)        / );
* demand - intermediate
loop((policy,policy_bau,r,s,g),  put  "id_pchg",        id_pchg.ts,        policy.tl,  policy_bau.tl,  r.tl,  "",      s.tl,  "",    "",      g.tl,  id_pchg(r,g,s,policy,policy_bau)      / );
loop((policy,policy_bau,g),      put  "id_pchg_g",      id_pchg_g.ts,      policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      g.tl,  id_pchg_g(g,policy,policy_bau)        / );
loop((policy,policy_bau,r),      put  "id_pchg_r",      id_pchg_r.ts,      policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    id_pchg_r(r,policy,policy_bau)        / );
loop((policy,policy_bau),        put  "id_pchg_us",     id_pchg_us.ts,     policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    id_pchg_us(policy,policy_bau)         / );
* demand - export and import
loop((policy,policy_bau,r,s),    put  "md_pchg",        md_pchg.ts,        policy.tl,  policy_bau.tl,  r.tl,  "",      s.tl,  "",    "",      "",    md_pchg(r,s,policy,policy_bau)        / );
loop((policy,policy_bau,r),      put  "md_pchg_r",      md_pchg_r.ts,      policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    md_pchg_r(r,policy,policy_bau)        / );
loop((policy,policy_bau),        put  "md_pchg_us",     md_pchg_us.ts,     policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    md_pchg_us(policy,policy_bau)         / );
loop((policy,policy_bau,r,s),    put  "mf_pchg",        mf_pchg.ts,        policy.tl,  policy_bau.tl,  r.tl,  "",      s.tl,  "",    "",      "",    mf_pchg(r,s,policy,policy_bau)        / );
loop((policy,policy_bau,r),      put  "mf_pchg_r",      mf_pchg_r.ts,      policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    mf_pchg_r(r,policy,policy_bau)        / );
loop((policy,policy_bau),        put  "mf_pchg_us",     mf_pchg_us.ts,     policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    mf_pchg_us(policy,policy_bau)         / );
* demand - consumption, investiment, and government
loop((policy,policy_bau,r,h),    put  "cd_pchg",        cd_pchg.ts,        policy.tl,  policy_bau.tl,  r.tl,  "",      "",    h.tl,  "",      "",    cd_pchg(r,h,policy,policy_bau)        / );
loop((policy,policy_bau,h),      put  "cd_pchg_h",      cd_pchg_h.ts,      policy.tl,  policy_bau.tl,  "",    "",      "",    h.tl,  "",      "",    cd_pchg_h(h,policy,policy_bau)        / );
loop((policy,policy_bau,r),      put  "cd_pchg_r",      cd_pchg_r.ts,      policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    cd_pchg_r(r,policy,policy_bau)        / );
loop((policy,policy_bau),        put  "cd_pchg_us",     cd_pchg_us.ts,     policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    cd_pchg_us(policy,policy_bau)         / );
loop((policy,policy_bau,r),      put  "i_pchg_r",       i_pchg_r.ts,       policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    i_pchg_r(r,policy,policy_bau)         / );
loop((policy,policy_bau),        put  "i_pchg_us",      i_pchg_us.ts,      policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    i_pchg_us(policy,policy_bau)          / );
loop((policy,policy_bau,r,pub),  put  "g_pchg",         g_pchg.ts,         policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    pub.tl,  "",    g_pchg(r,pub,policy,policy_bau)       / );
loop((policy,policy_bau,r),      put  "g_pchg_r",       g_pchg_r.ts,       policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    g_pchg_r(r,policy,policy_bau)         / );
loop((policy,policy_bau),        put  "g_pchg_us",      g_pchg_us.ts,      policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    g_pchg_us(policy,policy_bau)          / );
* tax revenue
loop((policy,policy_bau,r),      put  "trev_pchg_r",    trev_pchg_r.ts,    policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    trev_pchg_r(r,policy,policy_bau)      / );
loop((policy,policy_bau),        put  "trev_pchg_us",   trev_pchg_us.ts,   policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    trev_pchg_us(policy,policy_bau)       / );
* gross regional product
loop((policy,policy_bau,r),      put  "grpfc_pchg_r",   grpfc_pchg_r.ts,   policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    grpfc_pchg_r(r,policy,policy_bau)     / );
loop((policy,policy_bau),        put  "grpfc_pchg_us",  grpfc_pchg_us.ts,  policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    grpfc_pchg_us(policy,policy_bau)      / );
loop((policy,policy_bau,r),      put  "grpva_pchg_r",   grpva_pchg_r.ts,   policy.tl,  policy_bau.tl,  r.tl,  "",      "",    "",    "",      "",    grpva_pchg_r(r,policy,policy_bau)     / );
loop((policy,policy_bau),        put  "grpva_pchg_us",  grpva_pchg_us.ts,  policy.tl,  policy_bau.tl,  "",    "",      "",    "",    "",      "",    grpva_pchg_us(policy,policy_bau)      / );

putclose;

$label csv_save_cagr

$if %t_y%==2010 $goto wdrg_output

file results3 /'%save3_parameters%'/;
results3.pc = 5;
results3.nd = 6;
put results3;
* header
put                      "parameter",      "parameter.title", "simulation", "year", "region", "sector.aggregation", "sector", "household", "from", "value" /;
* compound annualized growth rate
* activity - production
loop((forloop,simtp,r,s),    put  "yd_cagr",         yd_cagr.ts,          forloop.tl,  simtp.tl,  r.tl,  "",      s.tl,  "",  "",  yd_cagr(r,s,simtp,forloop)        / );
loop((forloop,r,simtp),      put  "yd_cagr_r",       yd_cagr_r.ts,        forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",  "",  yd_cagr_r(r,simtp,forloop)        / );
loop((forloop,simtp),        put  "yd_cagr_us",      yd_cagr_us.ts,       forloop.tl,  simtp.tl,  "",    "",      "",    "",  "",  yd_cagr_us(simtp,forloop)         / );
loop((forloop,simtp,r,s),    put  "yxd_cagr",        yxd_cagr.ts,         forloop.tl,  simtp.tl,  r.tl,  "",      s.tl,  "",  "",  yxd_cagr(r,s,simtp,forloop)       / );
loop((forloop,r,simtp),      put  "yxd_cagr_r",      yxd_cagr_r.ts,       forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",  "",  yxd_cagr_r(r,simtp,forloop)       / );
loop((forloop,simtp),        put  "yxd_cagr_us",     yxd_cagr_us.ts,      forloop.tl,  simtp.tl,  "",    "",      "",    "",  "",  yxd_cagr_us(simtp,forloop)        / );
loop((forloop,simtp,r,s),    put  "yxf_cagr",        yxf_cagr.ts,         forloop.tl,  simtp.tl,  r.tl,  "",      s.tl,  "",  "",  yxf_cagr(r,s,simtp,forloop)       / );
loop((forloop,r,simtp),      put  "yxf_cagr_r",      yxf_cagr_r.ts,       forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",  "",  yxf_cagr_r(r,simtp,forloop)       / );
loop((forloop,simtp),        put  "yxf_cagr_us",     yxf_cagr_us.ts,      forloop.tl,  simtp.tl,  "",    "",      "",    "",  "",  yxf_cagr_us(simtp,forloop)        / );
loop((forloop,simtp,r,s),    put  "yt_cagr",         yt_cagr.ts,          forloop.tl,  simtp.tl,  r.tl,  "",      s.tl,  "",  "",  yt_cagr(r,s,simtp,forloop)        / );
loop((forloop,r,simtp),      put  "yt_cagr_r",       yt_cagr_r.ts,        forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",  "",  yt_cagr_r(r,simtp,forloop)        / );
loop((forloop,simtp),        put  "yt_cagr_us",      yt_cagr_us.ts,       forloop.tl,  simtp.tl,  "",    "",      "",    "",  "",  yt_cagr_us(simtp,forloop)         / );
* aggregate production
loop((forloop,simtp,r,s_a),  put  "yd_cagr_sa",      yd_cagr_sa.ts,       forloop.tl,  simtp.tl,  r.tl,  s_a.tl,  "",    "",  "",  yd_cagr_sa(r,s_a,simtp,forloop)   / );
loop((forloop,simtp,r,s_a),  put  "yxd_cagr_sa",     yxd_cagr_sa.ts,      forloop.tl,  simtp.tl,  r.tl,  s_a.tl,  "",    "",  "",  yxd_cagr_sa(r,s_a,simtp,forloop)  / );
loop((forloop,simtp,r,s_a),  put  "yxf_cagr_sa",     yxf_cagr_sa.ts,      forloop.tl,  simtp.tl,  r.tl,  s_a.tl,  "",    "",  "",  yxf_cagr_sa(r,s_a,simtp,forloop)  / );
loop((forloop,simtp,r,s_a),  put  "yt_cagr_sa",      yt_cagr_sa.ts,       forloop.tl,  simtp.tl,  r.tl,  s_a.tl,  "",    "",  "",  yt_cagr_sa(r,s_a,simtp,forloop)   / );
* gross regional product
loop((forloop,r,simtp),      put  "grpfc_cagr_r",    grpfc_cagr_r.ts,     forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",  "",  grpfc_cagr_r(r,simtp,forloop)     / );
loop((forloop,simtp),        put  "grpfc_cagr_us",   grpfc_cagr_us.ts,    forloop.tl,  simtp.tl,  "",    "",      "",    "",  "",  grpfc_cagr_us(simtp,forloop)      / );
loop((forloop,r,simtp),      put  "grpva_cagr_r",    grpva_cagr_r.ts,     forloop.tl,  simtp.tl,  r.tl,  "",      "",    "",  "",  grpva_cagr_r(r,simtp,forloop)     / );
loop((forloop,simtp),        put  "grpva_cagr_us",   grpva_cagr_us.ts,    forloop.tl,  simtp.tl,  "",    "",      "",    "",  "",  grpva_cagr_us(simtp,forloop)      / );
* total compound annualized growth rate
* activity - production
loop((forloop,r,s),          put  "yd_tcagr",        yd_tcagr.ts,         forloop.tl,  "",        r.tl,  "",      s.tl,  "",  "",  yd_tcagr(r,s,forloop)             / );
loop((forloop,r),            put  "yd_tcagr_r",      yd_tcagr_r.ts,       forloop.tl,  "",        r.tl,  "",      "",    "",  "",  yd_tcagr_r(r,forloop)             / );
loop((forloop),              put  "yd_tcagr_us",     yd_tcagr_us.ts,      forloop.tl,  "",        "",    "",      "",    "",  "",  yd_tcagr_us(forloop)              / );
loop((forloop,r,s),          put  "yxd_tcagr",       yxd_tcagr.ts,        forloop.tl,  "",        r.tl,  "",      s.tl,  "",  "",  yxd_tcagr(r,s,forloop)            / );
loop((forloop,r),            put  "yxd_tcagr_r",     yxd_tcagr_r.ts,      forloop.tl,  "",        r.tl,  "",      "",    "",  "",  yxd_tcagr_r(r,forloop)            / );
loop((forloop),              put  "yxd_tcagr_us",    yxd_tcagr_us.ts,     forloop.tl,  "",        "",    "",      "",    "",  "",  yxd_tcagr_us(forloop)             / );
loop((forloop,r,s),          put  "yxf_tcagr",       yxf_tcagr.ts,        forloop.tl,  "",        r.tl,  "",      s.tl,  "",  "",  yxf_tcagr(r,s,forloop)            / );
loop((forloop,r),            put  "yxf_tcagr_r",     yxf_tcagr_r.ts,      forloop.tl,  "",        r.tl,  "",      "",    "",  "",  yxf_tcagr_r(r,forloop)            / );
loop((forloop),              put  "yxf_tcagr_us",    yxf_tcagr_us.ts,     forloop.tl,  "",        "",    "",      "",    "",  "",  yxf_tcagr_us(forloop)             / );
loop((forloop,r,s),          put  "yt_tcagr",        yt_tcagr.ts,         forloop.tl,  "",        r.tl,  "",      s.tl,  "",  "",  yt_tcagr(r,s,forloop)             / );
loop((forloop,r),            put  "yt_tcagr_r",      yt_tcagr_r.ts,       forloop.tl,  "",        r.tl,  "",      "",    "",  "",  yt_tcagr_r(r,forloop)             / );
loop((forloop),              put  "yt_tcagr_us",     yt_tcagr_us.ts,      forloop.tl,  "",        "",    "",      "",    "",  "",  yt_tcagr_us(forloop)              / );
* aggregate production
loop((forloop,r,s_a),        put  "yd_tcagr_sa",     yd_tcagr_sa.ts,      forloop.tl,  "",        r.tl,  s_a.tl,  "",    "",  "",  yd_tcagr_sa(r,s_a,forloop)        / );
loop((forloop,r,s_a),        put  "yxd_tcagr_sa",    yxd_tcagr_sa.ts,     forloop.tl,  "",        r.tl,  s_a.tl,  "",    "",  "",  yxd_tcagr_sa(r,s_a,forloop)       / );
loop((forloop,r,s_a),        put  "yxf_tcagr_sa",    yxf_tcagr_sa.ts,     forloop.tl,  "",        r.tl,  s_a.tl,  "",    "",  "",  yxf_tcagr_sa(r,s_a,forloop)       / );
loop((forloop,r,s_a),        put  "yt_tcagr_sa",     yt_tcagr_sa.ts,      forloop.tl,  "",        r.tl,  s_a.tl,  "",    "",  "",  yt_tcagr_sa(r,s_a,forloop)        / );
* gross regional product
loop((forloop,r),            put  "grpfc_tcagr_r",   grpfc_tcagr_r.ts,    forloop.tl,  "",        r.tl,  "",      "",    "",  "",  grpfc_tcagr_r(r,forloop)          / );
loop((forloop),              put  "grpfc_tcagr_us",  grpfc_tcagr_us.ts,   forloop.tl,  "",        "",    "",      "",    "",  "",  grpfc_tcagr_us(forloop)           / );
loop((forloop,r),            put  "grpva_tcagr_r",   grpva_tcagr_r.ts,    forloop.tl,  "",        r.tl,  "",      "",    "",  "",  grpva_tcagr_r(r,forloop)          / );
loop((forloop),              put  "grpva_tcagr_us",  grpva_tcagr_us.ts,   forloop.tl,  "",        "",    "",      "",    "",  "",  grpva_tcagr_us(forloop)           / );

putclose;

$ontext
file results4 /'%save4_parameters%'/;
results4.pc = 5;
results4.nd = 6;
put results4;
put "parameter", "tp", "region", "sector", "household" , "inter.good", "value"/:
putclose;

file results5 /'%save5_parameters%'/;
results5.pc = 5;
results5.nd = 6;
put results5;
put "parameter", "tp", "region", "sector", "household" , "inter.good", "value"/:
putclose;

file results6 /'%save6_parameters%'/;
results6.pc = 5;
results6.nd = 6;
put results6;
put "parameter", "tp", "region", "sector", "household" , "inter.good", "value"/:
putclose;
$offtext

* ---------------------------------------------------------------------------- *
*        Coupled Model Output -- Files for the WDRG:
* ---------------------------------------------------------------------------- *

*               Label for WDRG

$label wdrg_output

$if not %s_wd%==1 $goto demand_response_output

*               Define region output names for WDRG

SET
         west_reg                DNDC regions (including ROUS)
                / Arizona, California, Colorado,
                  Idaho, Montana, Nevada,
                  NewMexico, Oregon, Utah,
                  Washington, Wyoming, ROUS /

        reg_map_ot(r,west_reg)      Cross-walk between DREM and DNDC regions (including ROUS)
                / arz . Arizona
                  cal . California
                  col . Colorado
                  ido . Idaho
                  mta . Montana
                  nmo . NewMexico
                  nva . Nevada
                  org . Oregon
                  uth . Utah
                  was . Washington
                  wyo . Wyoming
                  rus . ROUS /
;

*               Define WDRG specific parameters

PARAMETER
    prices        Crop rices for coupled WDRG Model - DNDC region r for DNDC crop sector s for time period tp under simulation sim
    production    Total crop production for coupled WDRG Model - DNDC region r for DNDC crop sector s for time period tp under simulation sim
;

prices(west_reg,crop_dndc,simtp,sim)     =  sum((reg_map_ot(r,west_reg), crop_map(crps,crop_dndc)), vyt_sim(r,crps,simtp, sim));
production(west_reg,crop_dndc,simtp,sim) =  sum((reg_map_ot(r,west_reg), crop_map(crps,crop_dndc)), yt_sim(r,crps,simtp,sim));

    DISPLAY prices, production;

execute_unload  '%save_wdrg_couple%'
                prices,
                production
;

execute_unload  '%save_wdrg_all%'
                vyt_sim, yt_sim,
                py_sim, pa_sim, pc_sim, pn_sim,
                yd_sim, yxd_sim, yxf_sim,
                yt_share_sa, yt_share_r,
                a_sim, md_sim, mf_sim,
                grpfc_sim_r, grpfc_sim_us,
                grpva_sim_r, grpva_sim_us,
                py_pchg, pa_pchg, pc_pchg, pn_pchg,
                yt_pchg, yd_pchg, yxd_pchg, yxf_pchg,
                a_pchg, md_pchg, mf_pchg,
                grpfc_pchg_r, grpfc_pchg_us,
                grpva_pchg_r, grpva_pchg_us
                yt_cagr, yd_cagr, yxd_cagr, yxf_cagr,
                grpfc_cagr_r, grpfc_cagr_us,
                grpva_cagr_r, grpva_cagr_us
;

* ---------------------------------------------------------------------------- *
*        Coupled Model Output -- Files for the Power System Model:
* ---------------------------------------------------------------------------- *

*               Label for demand response results

$label demand_response_output

$if not %s_dr%==1 $goto psm_output

execute_unload '%save_drm%'
                forloop, policy, policy_bau,
                r, zones, econr, econz_map, zones_map, econr_map, econzr_map,
                s, s_a, s_map,
                simtp, ts, tts,
                kstock_sim, lstock_sim, evok_sim, evol_sim, incadj_cal_sim,
                incadj_sim, rhscale_sim, ca_sim, hinv_sim, evolp_sim, evokp_sim,
                lpidx_sim,
                shk_dr_pchg_prod_sim, shk_dr_pchg_consum_sim,
                dre_prod_sim, dre_hh_sim,
                drate, irate, kprate, lprate, poprate,
                py_sim, pa_sim, pc_sim, pn_sim, pinv_sim,
                pgov_sim, pk_sim, pl_sim, pbt_sim,
                yd_sim, yd_sim_r, yd_sim_us,
                vyd_sim,
                yxd_sim, yxd_sim_r, yxd_sim_us,
                yxf_sim, yxf_sim_r, yxf_sim_us,
                yt_sim, yt_sim_r, yt_sim_us,
                yd_sim_sa, yxd_sim_sa, yxf_sim_sa, yt_sim_sa,
                kd_sim, kd_sim_r, kd_sim_us,
                ld_sim, ld_sim_r, ld_sim_us,
                btd_sim, btd_sim_r, btd_sim_us,
                id_sim, id_sim_g, id_sim_r, id_sim_us,
                trev_sim_r,
                cd_sim, cd_sim_r, cd_sim_us,
                a_sim,
                md_sim, md_sim_r, md_sim_us,
                mf_sim, mf_sim_r, mf_sim_us,
                cd_sim, cd_sim_h, cd_sim_r, cd_sim_us,
                g_sim, g_sim_r, g_sim_us,
                i_sim_r, i_sim_us,
                trev_sim_r, trev_sim_us,
                grpfc_sim_r, grpfc_sim_us,
                grpva_sim_r, grpva_sim_us,
                py_pchg, pa_pchg, pc_pchg, pn_pchg, pinv_pchg,
                pgov_pchg, pk_pchg, pl_pchg, pbt_pchg,
                yd_pchg, yd_pchg_r, yd_pchg_us,
                vyd_pchg,
                yxd_pchg, yxd_pchg_r, yxd_pchg_us,
                yxf_pchg, yxf_pchg_r, yxf_pchg_us,
                yt_pchg, yt_pchg_r, yt_pchg_us,
                yd_pchg_sa, yxd_pchg_sa, yxf_pchg_sa, yt_pchg_sa,
                a_pchg,
                kd_pchg, kd_pchg_r, kd_pchg_us,
                ld_pchg, ld_pchg_r, ld_pchg_us,
                btd_pchg, btd_pchg_r, btd_pchg_us,
                id_pchg, id_pchg_g, id_pchg_r, id_pchg_us,
                md_pchg, md_pchg_r, md_pchg_us,
                mf_pchg, mf_pchg_r, mf_pchg_us,
                cd_pchg, cd_pchg_h, cd_pchg_r, cd_pchg_us,
                i_pchg_r, i_pchg_us,
                g_pchg, g_pchg_r, g_pchg_us,
                trev_pchg_r, trev_pchg_us,
                grpfc_pchg_r, grpfc_pchg_us,
                grpva_pchg_r, grpva_pchg_us,
                yd_cagr,  yd_cagr_r,  yd_cagr_us,
                yxd_cagr, yxd_cagr_r, yxd_cagr_us,
                yxf_cagr, yxf_cagr_r, yxf_cagr_us,
                yt_cagr,  yt_cagr_r,  yt_cagr_us,
                yd_cagr_sa, yxd_cagr_sa, yxf_cagr_sa, yt_cagr_sa,
                grpfc_cagr_r, grpfc_cagr_us,
                grpva_cagr_r, grpva_cagr_us,
                yd_tcagr, yd_tcagr_r,  yd_tcagr_us,
                yxd_tcagr, yxd_tcagr_r, yxd_tcagr_us,
                yxf_tcagr, yxf_tcagr_r, yxf_tcagr_us,
                yt_tcagr, yt_tcagr_r, yt_tcagr_us,
                yd_tcagr_sa, yxd_tcagr_sa, yxf_tcagr_sa, yt_tcagr_sa,
                grpfc_tcagr_r, grpfc_tcagr_us,
                grpva_tcagr_r, grpva_tcagr_us
;

* ---------------------------------------------------------------------------- *
*        Coupled Model Output -- Files for the Power System Model:
* ---------------------------------------------------------------------------- *

*                Define PSM output label

$label psm_output

$if not %s_pw%==1 $goto end_script

*                Define PSM policy
SET
* Dynamic subsets
      psmpolicy(sim)
      psmbau(sim)

;

if ((switch_psm = 1 and
        switch_psmf = 0 and
        switch_dre = 0),

        psmpolicy(sim)$sim_psm(sim)    = YES;
        psmbau(sim)$sim_bau(sim)       = YES;

elseif (switch_psm = 1 and
        switch_psmf = 1 and
        switch_dre = 0),

        psmpolicy(sim)$sim_psmfull(sim) = YES;
        psmbau(sim)$sim_bau(sim)        = YES;

elseif (switch_psm = 1 and
         switch_psmf = 0 and
         switch_dre = 1),

       psmpolicy(sim)$sim_psmdrm(sim)  = YES;
       psmbau(sim)$sim_drm_sts(sim)    = YES;

elseif (switch_psm = 1 and
         switch_psmf = 1 and
         switch_dre = 1),

       psmpolicy(sim)$sim_psmfulldrm(sim)  = YES;
       psmbau(sim)$sim_drm_sts(sim)    = YES;

);

PARAMETER
        d_ely(r,sim,ssim) percentage change from the baseline - sectoral production for the domestic market - region r for sector electricity at time period 2050 under physical shock CNTF
        p_ely(r,sim,ssim) percentage change from the baseline - sectoral output price for the domestic market - region r for sector Electricity at time period 2050 under physical shock CNTF
        pa_ely(r,sim,ssim) percentage change from the baseline - armington aggregate price - region CA for sector Electricity at time period 2050 under physical shock CNTF

        d_ele_zones_pchg  percentage change from the baseline - sectoral production for the domestic market - zone zones for sector electricity at time period 2050 under physical shock CNTF
        d_ele_econr_pchg  percentage change from the baseline - sectoral production for the domestic market - state GenArea for sector electricity at time period 2050 under physical shock CNTF
        p_ele_zones_pchg  percentage change from the baseline - sectoral output price for the domestic market - zone zones for sector Electricity at time period 2050 under physical shock CNTF
        p_ele_econr_pchg  percentage change from the baseline - sectoral output price for the domestic market - state GenArea  for sector Electricity at time period 2050 under physical shock CNTF
        pa_ele_zones_pchg percentage change from the baseline - armington aggregate price - zone zones for sector Electricity at time period 2050 under physical shock CNTF
        pa_ele_econr_pchg percentage change from the baseline - armington aggregate price - state GenArea for sector Electricity at time period 2050 under physical shock CNTF

;

*SCALAR
*        d_elyCA      percentage change from the baseline - sectoral production for the domestic market - region r for sector electricity at time period 2050 under physical shock CNTF
*        d_elyROWECC  percentage change from the baseline - sectoral production for the domestic market - region r for sector electricity at time period 2050 under physical shock CNTF
*        p_elyCA      percentage change from the baseline - sectoral output price for the domestic market - region r for sector Electricity at time period 2050 under physical shock CNTF
*        p_elyROWECC  percentage change from the baseline - sectoral output price for the domestic market - region r for sector Electricity at time period 2050 under physical shock CNTF
*        pa_elyCA     percentage change from the baseline - armington aggregate price - region CA for sector Electricity at time period 2050 under physical shock CNTF
*        pa_elyROWECC percentage change from the baseline - armington aggregate price - region CA for sector Electricity at time period 2050 under physical shock CNTF

*;

LOOP(policy(sim),
    LOOP(policy_bau(ssim),
        d_ely(r,policy,policy_bau)$(yd_sim(r,'ele','%t_y%',policy_bau) <> 0)  = (yd_sim(r,'ele','%t_y%',policy)/yd_sim(r,'ele','%t_y%',policy_bau))-1;
        p_ely(r,policy,policy_bau)$(py_sim(r,'ele','%t_y%',policy_bau) <> 0)  = (py_sim(r,'ele','%t_y%',policy)/py_sim(r,'ele','%t_y%',policy_bau))-1;
        pa_ely(r,policy,policy_bau)$(pa_sim(r,'ele','%t_y%',policy_bau) <> 0) = (pa_sim(r,'ele','%t_y%',policy)/pa_sim(r,'ele','%t_y%',policy_bau))-1;
    );
);

LOOP(psmpolicy(sim),
    LOOP(psmbau(ssim),
*        d_elyCA      = (sum(r_map('CA',r),yd_sim(r,'ele','%t_y%',psmpolicy))/sum(r_map('CA',r),yd_sim(r,'ele','%t_y%',psmbau)))-1;
*        d_elyROWECC  = (sum(r_map('ROWECC',r),yd_sim(r,'ele','%t_y%',psmpolicy))/sum(r_map('ROWECC',r),yd_sim(r,'ele','%t_y%',psmbau)))-1;
*        p_elyCA      = (sum(r_map('CA',r),py_sim(r,'ele','%t_y%',psmpolicy))/sum(r_map('CA',r),py_sim(r,'ele','%t_y%',psmbau)))-1;
*        p_elyROWECC  = (sum(r_map('ROWECC',r),py_sim(r,'ele','%t_y%',psmpolicy))/sum(r_map('ROWECC',r),py_sim(r,'ele','%t_y%',psmbau)))-1;
*        pa_elyCA     = (sum(r_map('CA',r),pa_sim(r,'ele','%t_y%',psmpolicy))/sum(r_map('CA',r),pa_sim(r,'ele','%t_y%',psmbau)))-1;
*        pa_elyROWECC = (sum(r_map('ROWECC',r),pa_sim(r,'ele','%t_y%',psmpolicy))/sum(r_map('ROWECC',r),pa_sim(r,'ele','%t_y%',psmbau)))-1;

        d_ele_zones_pchg(zones)  = (sum(zones_map(zones,r),yd_sim(r,'ele','%t_y%',psmpolicy))/sum(zones_map(zones,r),yd_sim(r,'ele','%t_y%',psmbau)))-1;
        p_ele_zones_pchg(zones)  = (sum(zones_map(zones,r),py_sim(r,'ele','%t_y%',psmpolicy))/sum(zones_map(zones,r),py_sim(r,'ele','%t_y%',psmbau)))-1;
        pa_ele_zones_pchg(zones) = (sum(zones_map(zones,r),pa_sim(r,'ele','%t_y%',psmpolicy))/sum(zones_map(zones,r),pa_sim(r,'ele','%t_y%',psmbau)))-1;

        d_ele_econr_pchg(econr)  = (sum(econr_map(econr,r),yd_sim(r,'ele','%t_y%',psmpolicy))/sum(econr_map(econr,r),yd_sim(r,'ele','%t_y%',psmbau)))-1;
        p_ele_econr_pchg(econr)  = (sum(econr_map(econr,r),py_sim(r,'ele','%t_y%',psmpolicy))/sum(econr_map(econr,r),py_sim(r,'ele','%t_y%',psmbau)))-1;
        pa_ele_econr_pchg(econr) = (sum(econr_map(econr,r),pa_sim(r,'ele','%t_y%',psmpolicy))/sum(econr_map(econr,r),pa_sim(r,'ele','%t_y%',psmbau)))-1;

    );
);

DISPLAY
                psmpolicy, psmbau,
                d_ely, p_ely, pa_ely,
*                d_elyCA, p_elyCA, pa_elyCA,
*                d_elyROWECC, p_elyROWECC, pa_elyROWECC
                d_ele_zones_pchg, p_ele_zones_pchg, pa_ele_zones_pchg,
                d_ele_econr_pchg, p_ele_econr_pchg, pa_ele_econr_pchg
;

execute_unload '%save_psm_couple%'
                forloop, policy, psmpolicy, psmbau,
                r, zones, econr, econz_map, zones_map, econr_map, econzr_map,
                s, s_a, s_map,
                d_ely, p_ely, pa_ely,
*                d_elyCA, p_elyCA, pa_elyCA,
*                d_elyROWECC, p_elyROWECC, pa_elyROWECC
                d_ele_zones_pchg, p_ele_zones_pchg, pa_ele_zones_pchg,
                d_ele_econr_pchg, p_ele_econr_pchg, pa_ele_econr_pchg
;

execute_unload '%save_psm_all%'
                forloop, policy, policy_bau,
                psmpolicy, psmbau,
                r, zones, econr, econz_map, zones_map, econr_map, econzr_map,
                s, s_a, s_map,
                simtp, ts, tts,
                psmtfp_sim, psmtfp_fa_sim,
                kstock_sim, lstock_sim, evok_sim, evol_sim, incadj_cal_sim,
                incadj_sim, rhscale_sim, ca_sim, hinv_sim, evolp_sim, evokp_sim,
                lpidx_sim,
                psmtfp_fa_temp, psmtfp_temp,
                psmtfp_sim, psmtfp_fa_sim,
                shk_dr_pchg_prod_sim, shk_dr_pchg_consum_sim,
                drate, irate, kprate, lprate, poprate,
                py_sim, pa_sim, pc_sim, pn_sim, pinv_sim,
                pgov_sim, pk_sim, pl_sim, pbt_sim,
                yd_sim, yd_sim_r, yd_sim_us,
                vyd_sim,
                yxd_sim, yxd_sim_r, yxd_sim_us,
                yxf_sim, yxf_sim_r, yxf_sim_us,
                yt_sim, yt_sim_r, yt_sim_us,
                yd_sim_sa, yxd_sim_sa, yxf_sim_sa, yt_sim_sa,
                kd_sim, kd_sim_r, kd_sim_us,
                ld_sim, ld_sim_r, ld_sim_us,
                btd_sim, btd_sim_r, btd_sim_us,
                id_sim, id_sim_g, id_sim_r, id_sim_us,
                trev_sim_r,
                cd_sim, cd_sim_r, cd_sim_us,
                a_sim,
                md_sim, md_sim_r, md_sim_us,
                mf_sim, mf_sim_r, mf_sim_us,
                cd_sim, cd_sim_h, cd_sim_r, cd_sim_us,
                g_sim, g_sim_r, g_sim_us,
                i_sim_r, i_sim_us,
                trev_sim_r, trev_sim_us,
                grpfc_sim_r, grpfc_sim_us,
                grpva_sim_r, grpva_sim_us,
                py_pchg, pa_pchg, pc_pchg, pn_pchg, pinv_pchg,
                pgov_pchg, pk_pchg, pl_pchg, pbt_pchg,
                yd_pchg, yd_pchg_r, yd_pchg_us,
                vyd_pchg,
                yxd_pchg, yxd_pchg_r, yxd_pchg_us,
                yxf_pchg, yxf_pchg_r, yxf_pchg_us,
                yt_pchg, yt_pchg_r, yt_pchg_us,
                yd_pchg_sa, yxd_pchg_sa, yxf_pchg_sa, yt_pchg_sa,
                a_pchg,
                kd_pchg, kd_pchg_r, kd_pchg_us,
                ld_pchg, ld_pchg_r, ld_pchg_us,
                btd_pchg, btd_pchg_r, btd_pchg_us,
                id_pchg, id_pchg_g, id_pchg_r, id_pchg_us,
                md_pchg, md_pchg_r, md_pchg_us,
                mf_pchg, mf_pchg_r, mf_pchg_us,
                cd_pchg, cd_pchg_h, cd_pchg_r, cd_pchg_us,
                i_pchg_r, i_pchg_us,
                g_pchg, g_pchg_r, g_pchg_us,
                trev_pchg_r, trev_pchg_us,
                grpfc_pchg_r, grpfc_pchg_us,
                grpva_pchg_r, grpva_pchg_us,
                yd_cagr,  yd_cagr_r,  yd_cagr_us,
                yxd_cagr, yxd_cagr_r, yxd_cagr_us,
                yxf_cagr, yxf_cagr_r, yxf_cagr_us,
                yt_cagr,  yt_cagr_r,  yt_cagr_us,
                yd_cagr_sa, yxd_cagr_sa, yxf_cagr_sa, yt_cagr_sa,
                grpfc_cagr_r, grpfc_cagr_us,
                grpva_cagr_r, grpva_cagr_us,
                yd_tcagr, yd_tcagr_r,  yd_tcagr_us,
                yxd_tcagr, yxd_tcagr_r, yxd_tcagr_us,
                yxf_tcagr, yxf_tcagr_r, yxf_tcagr_us,
                yt_tcagr, yt_tcagr_r, yt_tcagr_us,
                yd_tcagr_sa, yxd_tcagr_sa, yxf_tcagr_sa, yt_tcagr_sa,
                grpfc_tcagr_r, grpfc_tcagr_us,
                grpva_tcagr_r, grpva_tcagr_us
;

* ----------------------------------------------------------------------------*
*       Exit program
* ----------------------------------------------------------------------------*

$label end_script

$exit

