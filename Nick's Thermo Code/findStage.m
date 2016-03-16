function findStage(t,vars,pressurefn,Vbottle,p_atm)
global stage;
switch stage
    case 1
        if vars(6) >= Vbottle
            stage=stage+1;
            vars(6)=Vbottle;
        end
    case 2
        [~,pstar]=pressurefn(t,vars);
        if pstar < p_atm
            stage=stage+1;
        end
    case 3
        p=pressurefn(t,vars);
        if p <= p_atm
            stage=stage+1;
        end
    case 4
        % This is the ballistic phase and does not end until the
        % integration ends.
end
end