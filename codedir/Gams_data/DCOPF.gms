$title  Power System Model - DCOPF model

* ------------------------------------------------------------------------------
*               select data based on target name
* ------------------------------------------------------------------------------

*       main result outputs

$set matout  "'results.gdx',  nse, toc_r, gen_r, coal_gen_r, ngas_gen_r, nse_r, dem_rd, nse, nse_id, shed_id, nse_ird, shed_ird, gen_gd, ele_export_id, ele_import_id, trans_zz2d";
$set matout2 "'results2.gdx', gen_ftype_rft, steam_gen_r, nuclear_gen_r, biogas_gen_r, geo_gen_r, petrol_gen_r, intertiel_gen_r, wind_gen_r, nse_rd";
$set matout3 "'results3.gdx', hydro_gen_r, pshydro_gen_r, motorload_gen_r, solar_gen_r, nse_r, nse, coal_gen_gr, ngas_gen_gr";
$set matout4 "'MarginalElectricityPrice.gdx', mc_id, mc_ird, mdem_id, mdem_ird, nse_id, nse_zd, gen_rd, tovc_coal_gr, tovc_ngas_gr, toc_coal_r, toc_ngas_r";
$set matout5 "'TransmissionImpacts.gdx',      PowerTransfer, MargCostPowerFlow, MargCostLineLimits,CNUS.L";

*        Mort's additional result outputs

*$set matout6 "'QuadFlows.gdx',demandq1,demandq2,demandq3,demandq4,genq1,genq2,genq3,genq4,netgenq1,netgenq2,netgenq3,netgenq4";
$set matout7 "'GenByBus.gdx',GenByBus";
$set matout8 "'OutageByBus.gdx',OutageByBus";

*       Joseph's additional results output

$set matout9  "'results4.gdx', gzones, geconr, Bustozone, ieconr, gi, gftype, gen_ftrd, fueluse_ftrd, toc_ftrd, tovc_ftr, mc_ird, mdem_ird, nse_ird, shed_ird";
$set matout10 "'results5.gdx', gen_grd, fueluse_grd, toc_grd, tovc_gr";
$set matout11 "'results6.gdx', mdem_ir, mdem_sh_ird, ele_export_ird, ele_import_ird, ele_flow_lijd";

* ------------------------------------------------------------------------------
*               set listing file options to limit it's size
* ------------------------------------------------------------------------------

$offsymxref
$offsymlist
$offdigit

option decimals = 8;

* ------------------------------------------------------------------------------
*               initialize sets and parameters.
* ------------------------------------------------------------------------------

$Include "TrialImportData.gms"

* ------------------------------------------------------------------------------
*               load unit commitment generator g on/off status
* ------------------------------------------------------------------------------

Parameter
        Genonoroff(g,d)         generator g turned on\off status from unit commitment model
;

$set gdxin UCtoDCOPF
$GDXIN %gdxin%
$LOAD Genonoroff
$GDXIN

* ------------------------------------------------------------------------------
*               define sets for DCOPF model
* ------------------------------------------------------------------------------

sets
        ii(line,i,i)            ac line l thermal limit and susceptance values
;

*       dynamic set for line l if x1

ii(line,i,j)$dtl(line,i,j,'x1') = Yes;

* ------------------------------------------------------------------------------
*               define scalars for DCOPF model
* ------------------------------------------------------------------------------

scalars
        CNS                      cost of non-served electricity / 1000 /
;

* ------------------------------------------------------------------------------
*               define additiona DCOPF model parameters
* ------------------------------------------------------------------------------

$ontext
*       keep in case we couple coal and natural gas generation parameters

a(g)=HeatRate(g)*FuelCost(g) + VOM(g);
a(g)$[(FuelType(g) = 1) and (GenArea(g)=2)]  = HeatRate(g) * FuelCost(g) + VOM(g);
a(g)$[(FuelType(g) = 5) and (GenArea(g)=2)]  = HeatRate(g) * FuelCost(g) + VOM(g);
a(g)$[(FuelType(g) = 6) and (GenArea(g)=2)]  = HeatRate(g) * FuelCost(g) + VOM(g);
a(g)$[(FuelType(g) = 1) and (GenArea(g)=1)]  = HeatRate(g) * FuelCost(g) + VOM(g);
a(g)$[(FuelType(g) = 5) and (GenArea(g)=1)]  = HeatRate(g) * FuelCost(g) + VOM(g);
a(g)$[(FuelType(g) = 6) and (GenArea(g)=1)]  = HeatRate(g) * FuelCost(g) + VOM(g);
$offtext

* ------------------------------------------------------------------------------
*               define coupled model parameters for DCOPF
* ------------------------------------------------------------------------------

scalars
        dem_scale              multipler to scale up or down the load depending on the model year - 2010 == 1  / 1 /
;

Parameter
        dem_mult(econr)        converted percentage change in demand multipler from economic regions to buses
        MultipliedDemand(i,d)  electricity demand after multiplying the demand multipler from the regional economic model
;

*       define demand change multipler for economic region r -- recall this is an annual percentage change from rem

dem_mult(econr) = 1 + DChange_DC(econr);

*       initiallize scaled demand at bus i for hour d

MultipliedDemand(i,d) = p1Demand(i,d)*dem_scale;

*       apply demand multipler to demand at bus i in economic region r for hour d

loop(ieconr(i,econr),

    MultipliedDemand(i,d) = p1Demand(i,d)*dem_scale*dem_mult(econr);
    
);

    display DChange_DC, dem_scale, dem_mult;

* ------------------------------------------------------------------------------
*               define DCOPF model variables
* ------------------------------------------------------------------------------

* can we combine the variables are in both to just positive variable
variable
        z                      total cost systen cost of electricity generation
        p(line,i,j,d)          electricity flow (power) on transmission line line at bus i to bus j for hour d
        Theta(i,d)             phase angle at bus i for hour d
;

positive variable
        x(g,d)                 electricity generation (MWH) of generator g for hour d
        CNUS(i,d)              non-served electricity (MWH) in zone zones for hour d
        xaboveminload(g,d)     electricity generation (MWH) required above minimum demand for stable load
        ShedGen(i,d)           shed generation (MWH) in zone zones for hour d
;

* ------------------------------------------------------------------------------
*               define DCOPF model variables
* ------------------------------------------------------------------------------

equations
        TC                     "total system cost of generation ($) - ob.j fun. - sum(startup cost, operation cost, and non-served electricity cost)"
        Demand(i,d)            total demand met (MWH) at bus i for hour demand hour d
        Powerflow1(line,i,j,d) electricity flow (power) (MWH) on transmission line line at bus i to bus j for hour d
        LB1(g,d)               relation between generation and minload above generation (MWH)
        LB2(g,d)               maximum generation above minimum load that is possible (MWH)
;

* ------------------------------------------------------------------------------
*         objective function equation
* ------------------------------------------------------------------------------

TC..           z =e= (sum((g,d),((a(g))*((x(g,d))))))+sum((i,d),((CNUS(i,d)+ShedGen(i,d))*(CNS)));

* ------------------------------------------------------------------------------
*         constraint equations
* ------------------------------------------------------------------------------

*       demand constraint

Demand(i,d)..  sum(g$ig(i,g),(x(g,d)))-sum((line,j)$ii(line,i,j),p(line,i,j,d))+
                      sum((line,j)$ii(line,j,i),p(line,j,i,d)) =e=((MultipliedDemand(i,d)))-CNUS(i,d) + ShedGen(i,d);

*       non-renewable generator constraints (FuelType == {1,5,6,7,14,16,17})

PowerFlow1(line,i,j,d)$(ii(line,i,j))..  p(line,i,j,d)         =e= [(Theta(i,d)-Theta(j,d))*1000*(dtl(line,i,j,'x1'))];
LB1(gggg,d)..                            x(gggg,d)             =e= ((Pmin(gggg)))*Genonoroff(gggg,d)+xaboveminload(gggg,d);
LB2(gggg,d)..                            xaboveminload(gggg,d) =l= ((Pmax(gggg))-Pmin(gggg))*Genonoroff(gggg,d);

* ------------------------------------------------------------------------------
*         constraint equations - fixed
* ------------------------------------------------------------------------------

*       fixing renewable generation availible for generator g in california for hour d

x.fx(g,d)$[[FuelType(g)=2]  and (GenArea(g)=1)] = ((Pmax(g)*AvailFinal(g,d)));
x.fx(g,d)$[[FuelType(g)=4]  and (GenArea(g)=1)] = ((Pmax(g)*AvailFinal(g,d)))*0.5;
x.fx(g,d)$[[FuelType(g)=8]  and (GenArea(g)=1)] = ((Pmax(g)*AvailFinal(g,d)));
x.fx(g,d)$[[FuelType(g)=9]  and (GenArea(g)=1)] = ((Pmax(g)*AvailFinal(g,d)));
x.fx(g,d)$[[FuelType(g)=15] and (GenArea(g)=1)] = ((Pmax(g)*AvailFinal(g,d)));

*       fixing renewable generation availible for generator g not in california for hour d
                           
x.fx(g,d)$[[FuelType(g)=2]  and not (GenArea(g)=1)] = ((Pmax(g)*AvailFinal(g,d)));
x.fx(g,d)$[[FuelType(g)=4]  and not (GenArea(g)=1)] = ((Pmax(g)*AvailFinal(g,d)));
x.fx(g,d)$[[FuelType(g)=8]  and not (GenArea(g)=1)] = ((Pmax(g)*AvailFinal(g,d)));
x.fx(g,d)$[[FuelType(g)=9]  and not (GenArea(g)=1)] = ((Pmax(g)*AvailFinal(g,d)));
x.fx(g,d)$[[FuelType(g)=15] and not (GenArea(g)=1)] = ((Pmax(g)*AvailFinal(g,d)));

*        fixing generation level for other generators

x.fx(g,d)$[FuelType(g)=3]  = 0;
x.fx(g,d)$[FuelType(g)=10] = 0;
x.fx(g,d)$[FuelType(g)=11] = 0;
x.fx(g,d)$[FuelType(g)=12] = Pmax(g);
x.fx(g,d)$[FuelType(g)=13] = 0;
x.fx(g,d)$[FuelType(g)=18] = 0;
x.fx(g,d)$[FuelType(g)=19] = 0;
x.fx(g,d)$[FuelType(g)=20] = 0;

*        fixing generation level from WBM (if GenLevel is greater than 0.01)
*        if the shock is enabled, generation == 0; else no effect

x.fx(gOff,d)$[GenLevel(gOff,d)>0.01] = 0;

*        fixing generation level for dynamic set for nuclear and coal retirement generators (set gg)
*        if the shock is enabled, generation == 0; else no effect

x.fx(gg,d) = 0;

* ------------------------------------------------------------------------------
*         constraint equations - lower and upper bounds
* ------------------------------------------------------------------------------

*Setting thermal line limits. Changed values for calibration.
*Possible update can be required if number of zones changes etc.

*        electricity flow for all lines

p.up(line,i,j,d)$ii(line,i,j) =  dtl(line,i,j,'lm')+300;
p.lo(line,i,j,d)$ii(line,i,j) = -dtl(line,i,j,'lm')-300;

*        electricity flow constrains for transmission lines connected to bus j in califorian (zones == 4)

p.up(line,i,j,d)$[ii(line,i,j) and Bustozone(j,'4')]  = p.up(line,i,j,d)+360;
p.lo(line,i,j,d)$[ii(line,i,j) and Bustozone(j,'4')]  = p.lo(line,i,j,d)-360;

*        electricity flow constrains for transmission lines connected to bus j in british columbia (zones == 10)

p.up(line,i,j,d)$[ii(line,i,j) and Bustozone(j,'10')] = p.up(line,i,j,d)+300;
p.lo(line,i,j,d)$[ii(line,i,j) and Bustozone(j,'10')] = p.lo(line,i,j,d)-300;

*        maximum and minimum values of the avraible Theta and initalizing the first theta to be zero
Theta.up(i,d)            =  0.5;
Theta.lo(i,d)            = -0.5;
Theta.fx(i,d)$[ORD(i)=1] =  0;

* ------------------------------------------------------------------------------
*               set CPLEX solver options
* ------------------------------------------------------------------------------

*       specifies the maximum time in seconds that a solver may run before it terminates

option ResLim=10000;

*       define model

Model DCOPF /all/;

*       solver call

Solve DCOPF using LP minimizing z;

* ------------------------------------------------------------------------------
*               display total generation
* ------------------------------------------------------------------------------

    display x.L;

* ------------------------------------------------------------------------------
*               define results.gdx output parameters
* ------------------------------------------------------------------------------

parameter
        toc_r(econr)                    total cost of electricity generation ($) - economic region r
        gen_r(econr)                    electricity generation (MWH) - economic region r
        nse_r(econr)                    non-served electricity (MWH) - economic region econr - added 5      
        coal_gen_r(econr)               coal electricity generation (MWH) - economic region r
        ngas_gen_r(econr)               "natural gas CT, CCGT, steam, turbine electricity generation (MWH) - economic region r"        
        gen_gd(g,d)                     electricity generation (MWH) - generator g - hour d
        dem_rd(econr,d)                 electricity demand (MWH) - economic region r - hour d  - added 5       
        nse                             non-served electricity (MWH) - added 5
        nse_id(i,d)                     non-served electricity (MWH) - bus i - hour d  - added 5
        shed_id(i,d)                    shed electricity (MWH) - bus i - hour d  - added 5
        nse_ird(i,econr,d)              unserved energy (MWh) - bus and hour - plus 5
        shed_ird(i,econr,d)             shed energy (MWh) - bus and hour - plus 5
        ele_export_id(i,d)              electricity exports (MWH) - bus i - hour d  - added 5
        ele_import_id(i,d)              electricity imports (MWH) - bus i - hour d  - added 5
        trans_zz2d(zones,zones2,d)      electricity transfer (MWH) - zones to zones2 - hour d
;

toc_r(econr)      = sum((g,d)$geconr(g,econr), x.l(g,d)*a(g));
gen_r(econr)      = sum((g,d)$geconr(g,econr), x.l(g,d));
nse_r(econr)      = sum((i,d)$[ieconr(i,econr)], CNUS.l(i,d)) + 5;
coal_gen_r(econr) = sum((g,d)$[geconr(g,econr) and  FuelType(g)=1], x.l(g,d));
ngas_gen_r(econr) = sum((g,d)$[geconr(g,econr) and (FuelType(g)=5 or
                                                    FuelType(g)=6 or
                                                    FuelType(g)=7)], x.l(g,d));                                                      
gen_gd(g,d)       = x.l(g,d);
dem_rd(econr,d)   = sum(i$ieconr(i,econr), MultipliedDemand(i,d)) + 5;                                                      
nse               = sum((i,d), CNUS.l(i,d)) + 5;
nse_id(i,d)       = CNUS.l(i,d) + 5;
shed_id(i,d)      = ShedGen.l(i,d) + 5;
nse_ird(i,econr,d)  = (CNUS.l(i,d) + 5)$[ieconr(i,econr)];
shed_ird(i,econr,d) = (ShedGen.l(i,d) + 5)$[ieconr(i,econr)] ;
ele_export_id(i,d)         = -sum((line,j)$ii(line,i,j),p.l(line,i,j,d)) + 5;
ele_import_id(i,d)         =  sum((line,j)$ii(line,j,i),p.l(line,j,i,d)) + 5;
trans_zz2d(zones,zones2,d) =  sum((line,i,j)$(linetozone(line,zones,zones2)),p.l(line,i,j,d));
 
    display toc_r, gen_r, dem_rd, coal_gen_r, ngas_gen_r, nse_r;

* ------------------------------------------------------------------------------
*               define results2.gdx output parameters
* ------------------------------------------------------------------------------

parameter
        gen_ftype_rft(econr,ft)         electricity generation (MWH) - fuel type ft - economic region econr 
        steam_gen_r(econr)              electricity generation (MWH) - steam turbine - economic region econr
        nuclear_gen_r(econr)            electricity generation (MWH) - nuclear - economic region econr
        biogas_gen_r(econr)             electricity generation (MWH) - bio-gass - economic region econr
        geo_gen_r(econr)                electricity generation (MWH) - geothermal - economic region econr
        petrol_gen_r(econr)             electricity generation (MWH) - petroleum - economic region econr
        intertiel_gen_r(econr)          electricity generation (MWH) - itertie - economic region econr
        wind_gen_r(econr)               electricity generation (MWH) - wind (onshore) - economic region econr
        nse_rd(econr,d)                 non-served electricity (MWH) - economic region econr - hour d  - added 5
;

gen_ftype_rft(econr,ft) = sum((g,d)$[geconr(g,econr) and  gftype(g,ft)], x.l(g,d));
steam_gen_r(econr)      = sum((g,d)$[geconr(g,econr) and  FuelType(g)= 7], x.l(g,d));
nuclear_gen_r(econr)    = sum((g,d)$[geconr(g,econr) and  FuelType(g)=12], x.l(g,d));
biogas_gen_r(econr)     = sum((g,d)$[geconr(g,econr) and  FuelType(g)=14], x.l(g,d));
geo_gen_r(econr)        = sum((g,d)$[geconr(g,econr) and  FuelType(g)=17], x.l(g,d));
petrol_gen_r(econr)     = sum((g,d)$[geconr(g,econr) and  FuelType(g)=16], x.l(g,d));
intertiel_gen_r(econr)  = sum((g,d)$[geconr(g,econr) and  FuelType(g)=10], x.l(g,d));
wind_gen_r(econr)       = sum((g,d)$[geconr(g,econr) and  FuelType(g)= 2], x.l(g,d));
nse_rd(econr,d)         = sum(i$[ieconr(i,econr)], CNUS.l(i,d)) + 5;

    display gen_ftype_rft, nse_rd;            

* ------------------------------------------------------------------------------
*               define results3.gdx output parameters
* ------------------------------------------------------------------------------

Parameter
        hydro_gen_r(econr)              electricity generation (MWH) - hydro - economic region econr
        pshydro_gen_r(econr)            electricity generation (MWH) - ps hydro - economic region econr
        solar_gen_r(econr)              electricity generation (MWH) - solar (all) - economic region econr
        motorload_gen_r(econr)          electricity generation (MWH) - motorload - economic region econr
        coal_gen_gr(g,econr)            electricity generation (MWH) - coal - economic region r
        ngas_gen_gr(g,econr)            electricity generation (MWH) - natural gas (CT and CCGT) - economic region r
;

hydro_gen_r(econr)     = sum((g,d)$[geconr(g,econr) and  FuelType(g)= 4], x.l(g,d));
pshydro_gen_r(econr)   = sum((g,d)$[geconr(g,econr) and  FuelType(g)=13], x.l(g,d));
motorload_gen_r(econr) = sum((g,d)$[geconr(g,econr) and  FuelType(g)=20], x.l(g,d));
solar_gen_r(econr)     = sum((g,d)$[geconr(g,econr) and (FuelType(g)= 8 or
                                                         FuelType(g)= 9 or
                                                         FuelType(g)=15)], x.l(g,d));                                                                  
coal_gen_gr(g,econr)   = sum(d$[geconr(g,econr)     and  FuelType(g) = 1], x.l(g,d));
ngas_gen_gr(g,econr)   = sum(d$[geconr(g,econr)     and (FuelType(g)= 5 or
                                                         FuelType(g)= 6)], x.l(g,d));

    display hydro_gen_r, pshydro_gen_r, solar_gen_r, motorload_gen_r,
            coal_gen_gr, ngas_gen_gr;

* ------------------------------------------------------------------------------
*               define MarginalElectricityPrice.gdx output parameters
* ------------------------------------------------------------------------------

Parameter
        gen_rd(econr,d)                 electricity generation (MWH) - economic region r - hour d   - added 5
        mc_id(i,d)                      marginal cost of electricity generation ($\MWH) - bus i - hour d  - added 5
        mc_ird(i,econr,d)               marginal cost of electricity generation ($\MWH) - bus i - economic region r - hour d  - added 5
        mdem_id(i,d)                    electricity demand (adjusted) (MWH) - bus i - hour d  - added 5
        mdem_ird(i,econr,d)             electricity demand (adjusted) (MWH) - bus i - economic region r - hour d  - added 5
        tovc_coal_gr(g,econr)           total variable cost for electricity generation ($\MWH) - coal in economic region r
        tovc_ngas_gr(g,econr)           total variable cost for electricity generation ($\MWH) - natural gas (CT and CCGT) in economic region r
        toc_coal_r(econr)               total cost of electricity generation ($) - coal - economic region r
        toc_ngas_r(econr)               total cost of electricity generation ($) - natural gas (CT: CCGT) - economic region r
        nse_zd(zones,d)                 non-served electricity (MWH) - psm zone zones - hour d  - added 5
;

gen_rd(econr,d)       = sum(g$[geconr(g,econr)], x.l(g,d));
mc_id(i,d)            =  Demand.m(i,d) + 5;
mc_ird(i,econr,d)     = (Demand.m(i,d) + 5)$[ieconr(i,econr)];
mdem_id(i,d)          =  MultipliedDemand(i,d) + 5;
mdem_ird(i,econr,d)   = (MultipliedDemand(i,d) + 5)$[ieconr(i,econr)];
tovc_coal_gr(g,econr) = a(g)$[geconr(g,econr) and  FuelType(g)=1];
tovc_ngas_gr(g,econr) = a(g)$[geconr(g,econr) and (FuelType(g)=5 or FuelType(g)=6)];
toc_coal_r(econr)     = sum((g,d)$[geconr(g,econr) and  FuelType(g)=1], x.l(g,d)*a(g));;
toc_ngas_r(econr)     = sum((g,d)$[geconr(g,econr) and (FuelType(g)=5 or
                                                        FuelType(g)=6)], x.l(g,d)*a(g));;
nse_zd(zones,d)       = sum(i$Bustozone(i,zones),CNUS.l(i,d))+5;      
            
    display gen_rd,       mc_ird,       mdem_ird,
            tovc_coal_gr, tovc_ngas_gr,
            toc_coal_r,   toc_ngas_r,   nse_zd;

* ------------------------------------------------------------------------------
*               define TransmissionImpacts.gdx output parameters
* ------------------------------------------------------------------------------

Parameter
        PowerTransfer(line,i,j,d)       electricity transfer (MWH) - line l at bus i and bus j - hour d - added 5
        MargCostPowerFlow(line,i,j,d)   electricity marginal cost of electricity transfer constraint ($\MWH) - line l at bus i and bus j - hour d
        MargCostLineLimits(line,i,j,d)  electricity marginal cost of electricity transfer limits ($\MWH) - line l at bus i and bus j - hour d
;

PowerTransfer(line,i,j,d)$ii(line,i,j) = p.l(line,i,j,d)+5;
MargCostPowerFlow(line,i,j,d)          = Powerflow1.M(line,i,j,d);
MargCostLineLimits(line,i,j,d)         = p.M(line,i,j,d);

    display PowerTransfer, MargCostPowerFlow, MargCostLineLimits;

* ------------------------------------------------------------------------------
*               define GenByBus.gdx output parameters
* ------------------------------------------------------------------------------

Parameters
        GenByBus(i,d)                     electricity generation - at bus i - hour d
        
;

GenByBus(i,d) = sum(g$ig(i,g),(x.L(g,d)));

* ------------------------------------------------------------------------------
*               define OutageByBus.gdx output parameters
* ------------------------------------------------------------------------------

Parameters
        OutageByBus(i,d)                 electricity outage - bus i - hour d
;

OutageByBus(i,d) = sum(gg$ig(i,gg),Pmax(gg)) + sum(gOff$[[ig(i,gOff)] and [GenLevel(gOff,d)>0.01] and not gg(gOff)],Pmax(gOff));

* ------------------------------------------------------------------------------
*               define results4.gdx output parameters
* ------------------------------------------------------------------------------

parameters
        gen_ftrd(ft,econr,d)            electricity generation (MWh) - fuel type ft - economic region r - hour d
        fueluse_ftrd(ft,econr,d)        fuel usage (HR*generation) (MMBTU) - fuel type ft - economic region r - hour d
        toc_ftrd(ft,econr,d)            total cost of generation($) - fuel type ft - economic region r - hour d
        tovc_ftr(ft,econr)              total variable cost of generation (HR*FC+VOM) ($\MWH) - fuel type ft - economic region r

;

gen_ftrd(ft,econr,d)     = sum(g$[geconr(g,econr) and gftype(g,ft)], x.l(g,d));
fueluse_ftrd(ft,econr,d) = sum(g$[geconr(g,econr) and gftype(g,ft)], x.l(g,d)*HeatRate(g));
toc_ftrd(ft,econr,d)     = sum(g$[geconr(g,econr) and gftype(g,ft)], x.l(g,d)*a(g));
tovc_ftr(ft,econr)       = sum(g$[geconr(g,econr) and gftype(g,ft)], a(g));

* ------------------------------------------------------------------------------
*               define results5.gdx output parameters
* ------------------------------------------------------------------------------

parameter
        gen_grd(g,econr,d)              electricity generation (MWh) - generator g - economic region r - hour d
        fueluse_grd(g,econr,d)          fuel usage (HR*generation) (MMBTU) - generator g - economic region r - hour d
        toc_grd(g,econr,d)              total cost of generation($) - generator g - economic region r - hour d
        tovc_gr(g,econr)                total variable cost of generation (HR*FC+VOM) ($\MWH) - generator g - economic region r
;

gen_grd(g,econr,d)     = (x.l(g,d))$geconr(g,econr);
fueluse_grd(g,econr,d) = (x.l(g,d)*HeatRate(g))$geconr(g,econr);
toc_grd(g,econr,d)     = (x.l(g,d)*a(g))$geconr(g,econr);
tovc_gr(g,econr)       = (a(g))$geconr(g,econr);

* ------------------------------------------------------------------------------
*               define results6.gdx output parameters
* ------------------------------------------------------------------------------

parameters
        mdem_ir(i,econr)                electiricty demand (MWh) - bus
        mdem_sh_ird(i,econr,d)          share of demand at hour d (MWh) - bus and hour
        ele_export_ird(i,econr,d)       electricity exports (MWH) - bus i - hour d  - added 5
        ele_import_ird(i,econr,d)       electricity imports (MWH) - bus i - hour d  - added 5
        ele_flow_lijd(line,i,j,d)       line energy transfer (MWh) - line bus bus and hour - plus 5
;

mdem_ir(i,econr)           = sum(d, MultipliedDemand(i,d))$[ieconr(i,econr)];
mdem_sh_ird(i,econr,d)$[mdem_ir(i,econr) gt 0] = (MultipliedDemand(i,d)$[ieconr(i,econr)] / mdem_ir(i,econr));
mdem_sh_ird(i,econr,d)$[mdem_ir(i,econr) eq 0] = 0;
ele_export_ird(i,econr,d)  = (-sum((line,j)$ii(line,i,j),p.l(line,i,j,d)) + 5)$[ieconr(i,econr)];
ele_import_ird(i,econr,d)  = ( sum((line,j)$ii(line,j,i),p.l(line,j,i,d)) + 5)$[ieconr(i,econr)];
ele_flow_lijd(line,i,j,d)$ii(line,i,j) =  p.l(line,i,j,d)+5;

* ------------------------------------------------------------------------------
*               export DCOPF statue for ParallelizingTrial.m
* ------------------------------------------------------------------------------

Scalar
        ModelStatus     model status for ParallelizingTrial.m: See GAMS Model and status code for more information.
        SolverStatus    model solver status for ParallelizingTrial.m: See GAMS Model and status code for more information.
;

ModelStatus  = DCOPF.modelstat;
SolverStatus = DCOPF.solvestat;

    display ModelStatus,SolverStatus;

execute_unload 'ModelStatusCodesOPF.gdx',ModelStatus,SolverStatus;

* ------------------------------------------------------------------------------
*               export all result parameters
* ------------------------------------------------------------------------------

execute_unload %matout%;
execute_unload %matout2%;
execute_unload %matout3%;
execute_unload %matout4%;
execute_unload %matout5%;
*execute_unload %matout6%;
execute_unload %matout7%;
execute_unload %matout8%;
execute_unload %matout9%;
execute_unload %matout10%;
execute_unload %matout11%;

* ------------------------------------------------------------------------------
*               extras
* ------------------------------------------------------------------------------

$ontext
Parameters demandq1(d), demandq2(d),  demandq3(d),   demandq4(d);
Parameters genq1(d), genq2(d),  genq3(d),   genq4(d);
Parameters netgenq1(d), netgenq2(d),  netgenq3(d),   netgenq4(d);

demandq1(d) = sum(i$[BusQuad(i) EQ 1], MultipliedDemand(i,d));
demandq2(d) = sum(i$[BusQuad(i) EQ 2], MultipliedDemand(i,d));
demandq3(d) = sum(i$[BusQuad(i) EQ 3], MultipliedDemand(i,d));
demandq4(d) = sum(i$[BusQuad(i) EQ 4], MultipliedDemand(i,d));

genq1(d) = sum(i$[BusQuad(i) EQ 1], sum(g$ig(i,g),(x.L(g,d)))   );
genq2(d) = sum(i$[BusQuad(i) EQ 2], sum(g$ig(i,g),(x.L(g,d)))   );
genq3(d) = sum(i$[BusQuad(i) EQ 3], sum(g$ig(i,g),(x.L(g,d))) );
genq4(d) = sum(i$[BusQuad(i) EQ 4], sum(g$ig(i,g),(x.L(g,d)))  );


netgenq1(d) = genq1(d) - demandq1(d);
netgenq2(d) = genq2(d) - demandq2(d);
netgenq3(d) = genq3(d) - demandq3(d);
netgenq4(d) = genq4(d) - demandq4(d);

Display demandq1,demandq2,demandq3,demandq4,genq1,genq2,genq3,genq4,netgenq1,netgenq2,netgenq3,netgenq4;
$offtext

* ------------------------------------------------------------------------------
*               define additional output parameters which are not exported
* ------------------------------------------------------------------------------
$ontext
* Kept this section for now
Parameter
        Generation(g,d)                      electricity generation of generator g for hour d
        CaliElectric                         electricity generation in CA
        RoweccElectric                       electricity generation in RO
*        CaliCoalShare
*        CaliGasShare
*        RoweccCoalShare
*        RoweccGasShare
        PriceCali                            electricity price in CA (TC\GEN)
        PriceRO                              electricity price in RO (TC\GEN)
        EnergyperBus(i)                      electricity demand at bus i
        DemandFraction(i,d)                  fraction of electicity demand at bus i for hour d
        MarginalCost(i,d)                    marginal cost of demand (from obj fun) at bus i for hour d
        WeightedAvgPriceEcon(zonesforEcon)   electricity demand weighted average annual price in economic regions
        UnservedEnergyCaliPerHourgams(d)     non-served electricity in CA for hour d
        UnservedEnergyROPerHourgams(d)       non-served electricity in RO for hour d
        DemandCaliperHourgams(d)             electricity demand at CA for hour d
        DemandROperHourgams(d)               electricity demand at RO for hour d
        UnservedEnergyFractionCaliperHour(d) share of non-served electricity of electiricy demand in CA for hour d  - added 5
        UnservedEnergyFractionROperHour(d)   share of non-served electricity of electiricy demand in RO for hour d  - added 5
        ShedEnergyValueSum                   shed electricity
        PowerTransfer2(line,i,j,d)           electricity transfer for line l at bus i and bus j for hour d - added 5
        ZoneTransfer(zones,d)                electricity export of zone zones for hour d


;

Generation(g,d) = x.l(g,d);

CaliElectric    = sum((g,d)$[GenArea(g)=1],(x.l(g,d)));
RoweccElectric  = sum((g,d)$[(GenArea(g)=2)],(x.l(g,d)));

PriceCali       = TotalCostCali/GenCali;
PriceRo         = TotalcostRo/GenRO;

EnergyperBus(i) = sum(d,MultipliedDemand(i,d));

DemandFraction(i,d)$[EnergyperBus(i) gt 0] = MultipliedDemand(i,d) / EnergyperBus(i);
DemandFraction(i,d)$[EnergyperBus(i) eq 0] = 0;

MarginalCost(i,d)         = Demand.m(i,d)+5;
WeightedAvgPriceEcon("1") = (sum((i,d)$ieconr(i,"1"),(DemandFraction(i,d)*Demand.m(i,d))))/44;
WeightedAvgPriceEcon("2") = (sum((i,d)$ieconr(i,"2"),(DemandFraction(i,d)*Demand.m(i,d))))/268;

UnservedEnergyCaliPerHourgams(d) = sum((i)$[ieconr(i,"1")],CNUS.l(i,d));
UnservedEnergyROPerHourgams(d)   = sum((i)$[ieconr(i,"2")],CNUS.l(i,d));

DemandCaliperHourgams(d) = sum((i)$[ieconr(i,"1")],MultipliedDemand(i,d));
DemandROperHourgams(d)   = sum((i)$[ieconr(i,"2")],MultipliedDemand(i,d));

UnservedEnergyFractionCaliperHour(d) = (UnservedEnergyCaliPerHourgams(d)/DemandCaliperHourgams(d))+5;
UnservedEnergyFractionROperHour(d)   = (UnservedEnergyROPerHourgams(d)/DemandROperHourgams(d))+5;

ShedEnergyValueSum = sum((i,d),ShedEnergyValue(i,d));

PowerTransfer2(line,i,j,d)$ii(line,i,j) = p.l(line,i,j,d)+5;

ZoneTransfer(zones,d) = -sum(i$Bustozone(i,zones),sum((line,j)$ii(line,i,j),p.l(line,i,j,d)))+sum(i$Bustozone(i,zones),sum((line,j)$ii(line,j,i),p.l(line,j,i,d)))+5;

    display CoalGenCaliPerGen,                 GasGenCaliPerGen,
            CoalGenROPerGen,                   GasGenROPerGen,
            PriceCali,                         PriceRO
            EnergyperBus,                      DemandFraction,
            MarginalCost,                      WeightedAvgPriceEcon
            DemandCaliperHour,                 DemandROperHour,
            UnservedEnergyFractionCaliperHour, UnservedEnergyFractionROperHour,
            ShedEnergyValueSum;
$offtext
