function [T,P] = compass_Tk(tIn,tParam)
      % This is T
      T = cell(size(tIn,1),1);
      % tParam.nd is the number of continuos
      % tParam.nx is the number of hidden state
      % tParam.dLinkMap, the link function between X and continuous variable
      for k=1:size(tIn,1)
        temp = ones(tParam.nc,size(tParam.xM,1));  
        for i=1:tParam.nc
          for j=1:size(tParam.xM,1)
              if tParam.cLinkMap(i,j)
                 temp(i,j)= tIn(k,tParam.cLinkMap(i,j));
              end
          end
        end
        T{k}=temp;
      end

      P = zeros(tParam.nc,tParam.nIn);  
      for i=1:tParam.nc
          for j=1:tParam.nIn
              if tParam.cConstantUpdate(i,j)
                 P(i,j)= 1;
              end
          end
      end
end