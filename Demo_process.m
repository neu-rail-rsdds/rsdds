
%%%%%You need to execute the following overall commands in sequence, 
%%%%%or view the results of each function block in combination with the paper.
%% GLRNNR  Generate initial two-dimensional saliency results at multiple scales
    GLRNNR;
%% Generate depth outlier detection results form depth map
    Dep_Outlier;
%% Recalculated 2D saliency results after adding depth saliency map
    GLRNNR_D;
%% Final nonlinear combination result
    Final_Fuse;
%% Evaluation of performance indicators for saliency results
    result_evaluate;