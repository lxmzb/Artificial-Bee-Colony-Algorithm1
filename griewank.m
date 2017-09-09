  function ObjVal = griewank(Chrom) %  Employed=Colony(1:(ABCOpts.ColonySize/2),:);50x5    ObjVal = griewank(Chrom,switch1);
   [Nind, Nvar] = size(Chrom);                     %[50,5]
   nummer = repmat(1:Nvar, [Nind 1]);
   ObjVal = sum(((Chrom.^2) / 4000)')' - prod(cos(Chrom ./ sqrt(nummer))')' + 1; %50ÐÐ1ÁÐ
  