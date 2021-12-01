%% Fitness function
function [Obj] = ProtocolOptimizerObjectiveCOHORTS_V2(x,PA,newpatients_C1)   

Obj = zeros(200,1);
parfor ii = 1:200%length(newpatients_C1)
    
    [solTreat,TotalDose] = SimulateComboDrugModelCOHORTS(PA,ii,newpatients_C1,x);

   % Objective function: Area under the tumour curve+area under treatment curve
    TimeSeries = linspace(0,PA.tf,1001); %Create 1000 evenly space points inside the time domain
    
    if solTreat.x(end)<PA.tf
        CumulativeTumourBurden = 10^10;
    else
        EvalSol = deval(solTreat,TimeSeries,1); %Evaluate the tumour solution at the collocation points
        CumulativeTumourBurden = trapz(TimeSeries,EvalSol); %Calculate the tumour AUC
    end
    %The objective function we are trying to minimize.
    Obj_vec(ii) = real(CumulativeTumourBurden) + TotalDose;

end

Obj = sum(Obj_vec);


end