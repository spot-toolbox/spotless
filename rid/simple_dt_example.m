load simple_dt_example;
data = ridata({Ytrain,Ytest},{Utrain,Utest},{Ttrain,Ttest});
[data,g] = states(data,'time-delay',0,0);

train = trials(data,1);
test = trials(data,2);

model = ripmodel(data,1,[7 1],g); % -- True system outside
                                  %    model class.

%model = lse_obj(model,data); %-- Equation Error.

dsel = unf_select(train,300);  %-- Representative subset of samples.

model = rie_obj(model,dsel);  %-- Local Robust Identification Error.

model = model.optimize();
val = model.validate(data);