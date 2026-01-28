function [QualityResults]=EnergyCycleGenerationTest(model,provenience,ATPMreaction)

%Gabriela Canto
%Jan 27th 2025



%Input:
%     - model = strcut array with S matrix
%     - provenience = bigg,carveme, modelseed
%     - ATPMreaction = ATP maintenance reaction

[~,ExCh]=findExcRxns(model);
model.lb(ExCh) = 0;
model.ub(ExCh) = 1000;

if (nargin>2)
    model.lb(findRxnIDs(model,ATPMreaction))=0;
    model.ub(findRxnIDs(model,ATPMreaction))=1000;
    
end

if contains(provenience,'bigg')
    mets={'atp';'ctp';'gtp';'utp';'itp';'nadh';'nadph';'fadh2';'fmnh2';'q8h2';'mql8';'2dmmql8';'accoa';'glu__L';'h'};
    cyto='_c';
    peri='_p';
    ex='_e';
    h2o_c='h2o_c';
    adp_c='adp_c';
    cdp_c='cdp_c';
    gdp_c='gdp_c';
    udp_c='udp_c';
    idp_c='idp_c';
    nad_c='nad_c';
    fad_c='fad_c';
    h_c='h_c';
    pi_c='pi_c';
    fmn_c='fmn_c';
    q8_c='q8_c';
    mqn8_c='mqn8_c';
    dmmq8_c='2dmmq8_c';
    ac_c='ac_c';
    coa_c='coa_c';
    akg_c='akg_c';
    nh4_c='nh4_c';
    h_p='h_p';
    nadp_c='nadp_c';
elseif contains(provenience,'carveme')
    mets={'atp';'ctp';'gtp';'utp';'itp';'nadh';'nadph';'fadh2';'fmnh2';'q8h2';'mql8';'2dmmql8';'accoa';'glu__L';'h'};
    cyto='[C_c]';
    peri='[C_p]';
    ex='[C_e]';
    
    h2o_c='h2o[C_c]';
    adp_c='adp[C_c]';
    cdp_c='cdp[C_c]';
    gdp_c='gdp[C_c]';
    udp_c='udp[C_c]';
    idp_c='idp[C_c]';
    h_c='h[C_c]';
    fad_c='fad[C_c]';
    pi_c='pi[C_c]';
    nad_c='nad[C_c]';
    fmn_c='fmn[C_c]';
    q8_c='q8[C_c]';
    mqn8_c='mqn8[C_c]';
    dmmq8_c='2dmmq8[C_c]';
    ac_c='ac[C_c]';
    coa_c='coa[C_c]';
    akg_c='akg[C_c]';
    nh4_c='nh4[C_c]';
    h_p='h[C_p]';
    nadp_c='nadp[C_c]';
elseif contains(provenience,'modelseed')
    mets={'cpd00002';'cpd00052';'cpd00038';'cpd00062';'cpd00068';'cpd00004';'cpd00005';'cpd00982';'cpd01270';'cpd15561';'cpd15499';'cpd15353';'cpd00022';'cpd00023';'cpd00067'};
    cyto='_c0[c0]';
    peri='_p0[p0]';
    ex='_e0[e0]';
    
    h2o_c='cpd00001_c0[c0]';
    adp_c='cpd00008_c0[c0]';
    cdp_c='cpd00096_c0[c0]';
    gdp_c='cpd00031_c0[c0]';
    udp_c='cpd00014_c0[c0]';
    idp_c='cpd00090_c0[c0]';
    h_c='cpd00067_c0[c0]';
    nad_c='cpd00003_c0[c0]';
    fad_c='cpd00015_c0[c0]';
    pi_c='cpd00009_c0[c0]';
    fmn_c='cpd00050_c0[c0]';
    q8_c='cpd15560_c0[c0]';
    mqn8_c='cpd15500_c0[c0]';
    dmmq8_c='cpd15352_c0[c0]';
    ac_c='cpd00029_c0[c0]';
    coa_c='cpd00010_c0[c0]';
    akg_c='cpd00024_c0[c0]';
    nh4_c='cpd00013_c0[c0]';
    h_p='cpd00067_c0[c0]';
    nadp_c='cpd00006_c0[c0]';
else
    print('Provenience not found')
end

model.c(:)=0;

for i=1:length(mets)
    model2=model;
    met=mets{i,1};
    met_ID=findMetIDs(model2,[met,cyto]);
    if met_ID==0
        QC_results{i,1}=[met,cyto];
        QC_results{i,2}='Met not in model';
        QC_results{i,3}='Met not in model';
        continue
    end
    if i==1
        
        %ATP
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={[met,cyto],h2o_c,adp_c,pi_c,h_c};
        stoichCoeffList=[-1 -1 1 1 1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
    
    elseif i==2
        %ctp
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={[met,cyto],h2o_c,cdp_c,pi_c,h_c};
        stoichCoeffList=[-1 -1 1 1 1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
       
    elseif i==3
        %gtp
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={[met,cyto],h2o_c,gdp_c,pi_c,h_c};
        stoichCoeffList=[-1 -1 1 1 1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
       
    elseif i==4
        %utp
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={[met,cyto],h2o_c,udp_c,pi_c,h_c};
        stoichCoeffList=[-1 -1 1 1 1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
       
    elseif i==5
        %itp
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={[met,cyto],h2o_c,idp_c,pi_c,h_c};
        stoichCoeffList=[-1 -1 1 1 1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
       
    elseif i==6
        %NADH
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={[met,cyto],h_c,nad_c};
        stoichCoeffList=[-1 1 1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
       
    elseif i==7
        %NADPH
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={[met,cyto],h_c,nadp_c};
        stoichCoeffList=[-1 1 1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
      
    elseif i==8
        %FADH2
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={[met,cyto],h_c,fad_c};
        stoichCoeffList=[-1 2 1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
      
    elseif i==9
        %FMNH2
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={[met,cyto],h_c,fmn_c};
        stoichCoeffList=[-1 2 1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
     
    elseif i==10
        %Q8H2
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={[met,cyto],h_c,q8_c};
        stoichCoeffList=[-1 2 1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
    
    elseif i==11
        %MQL8
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={[met,cyto],h_c,mqn8_c};
        stoichCoeffList=[-1 2 1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
    
    elseif i==12
        %DMMQL8
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={[met,cyto],h_c,dmmq8_c};
        stoichCoeffList=[-1 2 1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
     
    elseif i==13
        %ACCOA
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={[met,cyto],h2o_c,h_c,ac_c,coa_c};
        stoichCoeffList=[-1 -1 1 1 1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
      
    elseif i==14
        %GLU
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={[met,cyto],h2o_c,akg_c,nh4_c,h_c};
        stoichCoeffList=[-1 -1 1 1 2];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
      
    elseif i==15
        %PROTON
        rxnName={[met,'QC'],[met,'QC']};
        metaboliteList={h_p,[met,cyto]};
        stoichCoeffList=[-1 1];
        revFlag=false;
        lowerBound=0;
        upperBound=1000;
        objCoeff=0;
        subSystem='QC analysis';
        grRule='';
        geneNameList='';
        systNameList='';
        checkDuplicate=true;
        confScores={1};
        rxnECNumbers='';
        rxnKEGGID='';
        rxnReferences='';
        rxnNotes='';
        model2.csense(length(model2.mets),1)='E';
       
    end
    model2=addReaction_EGC(model2,rxnName,metaboliteList,stoichCoeffList,revFlag,lowerBound,upperBound,objCoeff,subSystem,grRule,geneNameList,systNameList,checkDuplicate,confScores);  
    FBA=optimizeCbModel(model2);
    
    if FBA.f >0
       QC_results{i,1}=met;
       QC_results{i,2}=FBA.f;
       QC_results{i,3}=1;
    else
       QC_results{i,1}=met;
       QC_results{i,2}=FBA.f;
       QC_results{i,3}=0;
       
    end
    QualityResults=cell2table(QC_results,"VariableNames",["Met","Flux","EnergyProduction"]);
    
end



end