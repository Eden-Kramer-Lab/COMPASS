function []= compass_run_models(variables)
%% Store compass_em output
for i = 1:length(variables) %for as many variable sets given, a new field containing all compass_em output parameters will be created
    this = char(variables(i));
    load(this)
    [rXSmt,rSSmt,Param,rXPos,rSPos,ML,EYn,EYb,rYn,rYb] = compass_em(DISTR,Uk,In,Ib,Yn,Yb,Param,obs_valid);
    answer = {'rXSmt','rSSmt','Param','rXPos','rSPos','ML','EYn','EYb','rYn','rYb';rXSmt,rSSmt,Param,rXPos,rSPos,ML,EYn,EYb,rYn,rYb};
    save(['model_result' num2str(i)],'answer')
end
end
