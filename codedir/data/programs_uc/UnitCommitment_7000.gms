$title  Power System Model - Unit Commitment Program with a 7000 limit

* ------------------------------------------------------------------------------
$ontext
We are building a Unit Commitment model which will help us understand the 
generators to commit and use that information in the DC power flow equations.
$offtext
* ------------------------------------------------------------------------------

* ------------------------------------------------------------------------------
* 		load model data
* ------------------------------------------------------------------------------  

$Include "TrialImportData.gms"

* ----------------------------------------------------------------------------*
*       set system options:
* ----------------------------------------------------------------------------*

*       number of decimals

option decimals = 8;

* ------------------------------------------------------------------------------
* 		define coupled model parameters for the unit commitment model 
* ------------------------------------------------------------------------------ 


parameter
*   parameters from the economic model
*        DemandChangeCali     coupled Model Parameter (percentage change in demand -annually- in CA)
*        DemandChangeRO       coupled Model Parameter (percentage change in demand -annually- in ROWECC)
    
*   paramters translasted from REM to PSM   
        DemandMult(zones)    converted percentage change in demand multipler from economic regions to zones
;

*DemandChangeCali = DChange('1'); 
*DemandChangeRO   = DChange('2');

*DemandMult(zones)$ZoneLocation(zones,'1') = 1 + DemandChangeCali;
*DemandMult(zones)$ZoneLocation(zones,'2') = 1 + DemandChangeRO;

DemandMult(zones) = 1;
DemandMult(zones) = 1 + DChange(zones);

    display DChange, DemandMult;

* ------------------------------------------------------------------------------
* 		define intial parameters for the unit commitment model 
* ------------------------------------------------------------------------------ 

parameter
*   model parameters    
        Spinreserves(d)      spinning reserves at hour d as 1%

;

Spinreserves(d) = sum(zones, LoadDemand(zones,d))*0.01;

$ontext
*       keep in case we couple coal and natural gas generation parameters
a(g)                                         = HeatRate(g) * FuelCost(g) + VOM(g);
a(g)$[(FuelType(g) = 1) and (GenArea(g)=2)]  = HeatRate(g) * FuelCost(g) + VOM(g);
a(g)$[(FuelType(g) = 5) and (GenArea(g)=2)]  = HeatRate(g) * FuelCost(g) + VOM(g);
a(g)$[(FuelType(g) = 6) and (GenArea(g)=2)]  = HeatRate(g) * FuelCost(g) + VOM(g);
a(g)$[(FuelType(g) = 1) and (GenArea(g)=1)]  = HeatRate(g) * FuelCost(g) + VOM(g);
a(g)$[(FuelType(g) = 5) and (GenArea(g)=1)]  = HeatRate(g) * FuelCost(g) + VOM(g);
a(g)$[(FuelType(g) = 6) and (GenArea(g)=1)]  = HeatRate(g) * FuelCost(g) + VOM(g);
$offtext

* ------------------------------------------------------------------------------
* 		define variables for the Unit Commitment Program
* ------------------------------------------------------------------------------

variables
        z                                    total cost systen cost of electricity generation
;                                            
                                             
positive variables                           
        xTurnoff(g,d)                        binary variable (0\1) - determine generator g off status for hour d --- move down to binary variable?
        CNS(zones,d)                         non-served electricity (MWH) in zone zones for hour d
        x(g,d)                               electricity generation (MWH) of generator g for hour d
        xaboveminload(g,d)                   electricity generation (MWH) required above minimum demand for stable load
        Lineflow(lines,d)                    electricity flow (MWH) over zone transmission line lines for hour d
        ShedGen(zones,d)                     shed generation (MWH) in zone zones for hour d
;

binary variables
        xTurnon(g,d)                         binary variable (0\1) - determine generator g on status for hour d
        xonoroff(g,d)                        binary variable (0\1)- indicator of the generator g  up\down status for hour d

;

* ------------------------------------------------------------------------------
* 		define alias
* ------------------------------------------------------------------------------

*       hour alias

alias(d,dd);

* ------------------------------------------------------------------------------
* 		define Scalars for the Unit Commitment Program
* ------------------------------------------------------------------------------ 

Scalar
        CNUS                                non-served electricity cost ($) /1000/
;

* ------------------------------------------------------------------------------
* 		define Equations for the Unit Commitment Program
* ------------------------------------------------------------------------------

Equations
        TotalCost                           total system cost of generation ($) - ob.j fun. - sum(startup cost: operation cost: and non-served electricity cost)
        Demand(zones,d)                     total demand met in zone zones for hour demand hour d
        PoweraboveMinGen(g,d)               electricity generation required above the minimum load capacity
        RelBetweenMinGenandaboveMinGen(g,d) relationship equation between MinGen and AboveminGen
        PlantStatus(g,d)                    on\off status of generator g for hour d
        rampconstraintup(g,d)               ramping up constraint of generator g for hour d
        rampconstraintdown(g,d)             ramping down constraint of generator g for hour d
        enforceminuptime(g,d)               uptime of the generator g for hour d
        enforcemindowntime(g,d)             downtime of generator g for hour d
        meetreserves(d)                     minimum system reserve requirement
;

* ------------------------------------------------------------------------------
*	  objective function equation
* ------------------------------------------------------------------------------

TotalCost..     z=e=sum((zones,d),((CNS(zones,d)+ShedGen(zones,d))*CNUS))+sum((g,d),(x(g,d)*a(g)))+ sum((g,d),(Startupcost(g)*xTurnon(g,d)));

* ------------------------------------------------------------------------------
*	  constraint equations
* ------------------------------------------------------------------------------

*       demand constraint
Demand(zones,d)..                        sum(g$Gentozone(g,zones),(x(g,d))) - sum(lines$LineSources(lines,zones), Lineflow(lines,d))+
                                                    sum(lines$LineSink(lines,zones),Lineflow(lines,d)) =e= (LoadDemand(zones,d)*DemandMult(zones))-Cns(zones,d) +ShedGen(zones,d);

*       non-renewable generator constraints (FuelType == {1,5,6,7,14,16,17})
PoweraboveMinGen(gggg,d)..               xaboveminload(gggg,d) =l= (((Pmax(gggg))-Pmin(gggg)))*xonoroff(gggg,d);
RelBetweenMinGenandaboveMinGen(gggg,d).. x(gggg,d)             =e= xaboveminload(gggg,d)+((Pmin(gggg)*xonoroff(gggg,d)));
PlantStatus(gggg,d)..                    xonoroff(gggg,d)      =e= xonoroff(gggg,d-1)$[Ord(d)>1]+xTurnon(gggg,d)-xTurnoff(gggg,d);
rampconstraintup(gggg,d)$[ORD(d)>1]..    xaboveminload(gggg,d)-xaboveminload(gggg,d-1) =l= Ramp(gggg);
rampconstraintdown(gggg,d)$[ORD(d)>1]..  xaboveminload(gggg,d-1)-xaboveminload(gggg,d) =l= Ramp(gggg) ;
enforceminuptime(gggg,d) ..              xonoroff(gggg,d)      =g= sum(dd$[ORD(dd)<=ORD(d) and ORD(dd)>(ORD(d)-Uptime(gggg))],xTurnon(gggg,dd));
enforcemindowntime(gggg,d) ..            1-xonoroff(gggg,d)    =g= sum(dd$[ORD(dd)<=ORD(d) and ORD(dd)>(ORD(d)-Downtime(gggg))],xTurnoff(gggg,dd));

*       spinning reserve
meetreserves(d)..                        Spinreserves(d) =l= sum(gggg,xonoroff(gggg,d)*Pmax(gggg)-x(gggg,d));

* ------------------------------------------------------------------------------
*	  constraint equations - lower and upper bounds
* ------------------------------------------------------------------------------

*        zone transimission line constrain (maximum line flow between the zones)

Lineflow.up(lines,d) = LineCap(lines);


*        maximum generation for non-renewable generator g

x.up(gggg,d)         = Pmax(gggg);

* ------------------------------------------------------------------------------
*	  constraint equations - fixed
* ------------------------------------------------------------------------------

*       renewable energy availible for generator g

x.fx(g,d)$[FuelType(g)=2]  = ((Pmax(g)*AvailFinal(g,d)));
x.fx(g,d)$[FuelType(g)=3]  = 0;
*x.fx(g,d)$[[FuelType(g)=4] and (GenArea(g)=2)] = ((Pmax(g)*AvailFinal(g,d)))*0.5;
*x.fx(g,d)$[[FuelType(g)=4] and (GenArea(g)=1)] = ((Pmax(g)*AvailFinal(g,d)))*0.5;
x.fx(g,d)$[FuelType(g)=4]  = ((Pmax(g)*AvailFinal(g,d)))*0.5;
x.fx(g,d)$[FuelType(g)=8]  = ((Pmax(g)*AvailFinal(g,d)));
x.fx(g,d)$[FuelType(g)=9]  = ((Pmax(g)*AvailFinal(g,d)));
x.fx(g,d)$[FuelType(g)=10] = 0;
x.fx(g,d)$[FuelType(g)=11] = 0;
x.fx(g,d)$[FuelType(g)=12] = Pmax(g);
x.fx(g,d)$[FuelType(g)=13] = 0;
x.fx(g,d)$[FuelType(g)=15] = ((Pmax(g)*AvailFinal(g,d)));
x.fx(g,d)$[FuelType(g)=18] = 0;
x.fx(g,d)$[FuelType(g)=19] = 0;
x.fx(g,d)$[FuelType(g)=20] = 0;

*        fixing on\off status of renewable and nuclear plants.

xonoroff.fx(g,d)$[FuelType(g)=2]  = 1;
xonoroff.fx(g,d)$[FuelType(g)=3]  = 0;
xonoroff.fx(g,d)$[FuelType(g)=4]  = 1;
xonoroff.fx(g,d)$[FuelType(g)=8]  = 1;
xonoroff.fx(g,d)$[FuelType(g)=9]  = 1;
xonoroff.fx(g,d)$[FuelType(g)=10] = 0;
xonoroff.fx(g,d)$[FuelType(g)=11] = 0;
xonoroff.fx(g,d)$[FuelType(g)=12] = 1;
xonoroff.fx(g,d)$[FuelType(g)=13] = 0;
xonoroff.fx(g,d)$[FuelType(g)=15] = 1;
xonoroff.fx(g,d)$([(FuelType(g)=18) or (FuelType(g)=19) or (FuelType(g)=20) ]) = 0;

*        off status generators from WBM, if GenLevel is greater than 0.01

x.fx(gOff,d)$[GenLevel(gOff,d)>0.01]        = 0;
xonoroff.fx(gOff,d)$[GenLevel(gOff,d)>0.01] = 0;

*        off status generators in the dynamic set for nuclear and coal retirement generators (set gg)

xonoroff.fx(gg,d) = 0;
x.fx(gg,d)        = 0;

* ------------------------------------------------------------------------------
* 		set CPLEX solver options
* ------------------------------------------------------------------------------

*       specifies the maximum time in seconds that a solver may run before it terminates
      
option RESLIM = 7000;

*       specifies a relative termination tolerance

option optcr = 0.000001;

*       define model

Model UnitCommitment /all/;

*       additional model and solve statements

*UnitCommitment.holdfixed = 1 ;
UnitCommitment.dictfile=0;
UnitCommitment.optfile=1;
UnitCommitment.scaleopt=1;

UnitCommitment.holdFixed=1;
UnitCommitment.savepoint = 2;

*        export model results to speed up next unit commitment model run

execute_loadpoint 'UnitCommitment_p1';

*       solver call

Solve UnitCommitment using mip minimizing z;

* ------------------------------------------------------------------------------
* 		define unit commitment model outputs
* ------------------------------------------------------------------------------

scalar
        ModelStatus       model status for ParallelizingTrial.m: See GAMS Model and status code for more information.
        SolverStatus      model solver status for ParallelizingTrial.m: See GAMS Model and status code for more information.
;

parameter
        pMIPGAP          MIP Gap value (different between estimated objective and actual objective value) 
        pGAPPERC         percentage change of objective estimate from objective value
        Genonoroff(g,d)  parameter to tell DCOPF which generators are turned on and off
;


*       output for Matlab

ModelStatus      = UnitCommitment.modelstat;
SolverStatus     = UnitCommitment.solvestat;
pMIPGAP          = UnitCommitment.objEst - UnitCommitment.objVal;
pGAPPERC         = pMIPGAP / UnitCommitment.objVal;

*       output for DCOPF model

Genonoroff(g,d) = xonoroff.l(g,d);

    display ModelStatus,SolverStatus, pMIPGAP, pGAPPERC;

* ------------------------------------------------------------------------------
* 		export unit commitment model outputs
* ------------------------------------------------------------------------------

*       output for ParallelizingTrial.m

execute_unload 'ModelStatusCodes.gdx',ModelStatus,SolverStatus,pMIPGAP,pGAPPERC;

*       output for DCOPF

execute_unload 'UCtoDCOPF.gdx',Genonoroff ;
