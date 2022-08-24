function [res, est, optOut] = compare_methods(SRI, HSI, MSI, DegMat, ds, methods)
  res = cell(size(methods,1), 6);
  est = cell(size(methods,1),1);
  optOut = cell(size(methods,1),1);
  
  for i=1:size(methods,1)
    
    % Initialize parameters and find name
    opt = struct();
    blockstr = '';
    if size(methods,2) > 3 && ~isempty(methods{i,4})
        if strcmp(methods{i,1},'HySure')||strcmp(methods{i,1},'CNMF')||strcmp(methods{i,1},'FuVar') 
            opt.P = methods{i,4};
            blockstr = sprintf('P = %d', opt.P);
        elseif strcmp(methods{i,1},'CT-STAR') 
            opt.R2 = methods{i,4};
            blockstr = sprintf('R = ');
        elseif strcmp(methods{i,1},'CB-STAR') 
            opt.R2 = methods{i,4}{1};
            opt.initOpt = methods{i,4}{2};
            blockstr = sprintf('R = ');
        elseif strcmp(methods{i,1},'LTMR') 
            opt.p = methods{i,4}{1};
            opt.lambda = methods{i,4}{2};
            blockstr = sprintf('R = ');
        end
    end
    if strcmp(methods{i,1},'HySure') 
        methname = sprintf('%s %s%s', methods{i,1}, blockstr);
    else
        methname = sprintf('%s %s%s', methods{i,1}, blockstr, methods{i,3});
    end
    fprintf('Running method %s\n', methname);
    
    Params = sort(fieldnames(DegMat));  
    ParamStr = sprintf('DegMat.%s, ', Params{:});
    tic,
    % eval(sprintf('[Y_hat,~] = %s(HSI, MSI, %s %s, opt);', ...
    %         methods{i,2}, ParamStr,  methods{i,3}));
    eval(sprintf('[Y_hat,optOut{i}] = %s(HSI, MSI, %s %s, opt);', ...
            methods{i,2}, ParamStr,  methods{i,3}));
    time = toc;
    
    res(i,:) = [methname, compute_metrics(Y_hat,SRI,ds), time];
    est{i} = Y_hat;
  end       
end



