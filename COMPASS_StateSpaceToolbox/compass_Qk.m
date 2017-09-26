function [Q,P] = compass_Qk(tIn,tParam)
      % This is Qt
      Q = cell(size(tIn,1),1);
      P = cell(size(tIn,1),1);
      % tParam.nd is the number of continuos
      % tParam.nx is the number of hidden state
      % tParam.dLinkMap, the link function between X and continuous variable
      for k=1:size(tIn,1)
        temp = ones(tParam.nd,size(tParam.xM,1));  
        for i=1:tParam.nd
          for j=1:size(tParam.xM,1)
              if tParam.dLinkMap(i,j)
                 temp(i,j)= tIn(k,tParam.dLinkMap(i,j));
              end
          end
        end
        Q{k}=temp;
      end
      P = zeros(tParam.nd,tParam.nIb);  
      for i=1:tParam.nd
          for j=1:tParam.nIb
              if tParam.dConstantUpdate(i,j)==1
                 P(i,j)= 1;
              end
          end
      end
      
end