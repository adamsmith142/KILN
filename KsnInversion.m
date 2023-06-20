function inversion = KsnInversion(S, DEM, A, grid_size, alpha, mn)

A0 = 1e6;

%x y and z of all the stream nodes in order from highest elevation to
%lowest
x = S.x(S.ixc);
y = S.y(S.ixc);
z = DEM.Z(S.IXgrid(S.ixc));
z = double(z);

%calculate chi values
C = chitransform(S, A, 'mn', mn);
chi = C(S.ixc);

%make base level chi value 0
chi(length(chi)) = 0;

%produce a linear list of IDs
give = zeros([length(S.ix), 1]);
for i = 1:length(S.ix)
    give(i) = i;
end

%produce list of receiving IDs for the linear list
rec = zeros([length(S.ix), 1]);
for i = 1:length(S.ix)
    id = S.ixc(i);
    k = find(S.ix == id);
    if isempty(k)
        k = 0;
        rec(i) = k;
    else
        rec(i) = k;
   
        end 
end


%baselevel correction
z0 = zerobaselevel(S, DEM);
z0 = z0(S.ixc);

%calculate distance between nodes
d = distance(S, 'node_to_node');
d = d(S.ix);

%calculate upstream area
area = A.Z(S.IXgrid(S.ix))*(A.cellsize^2);

%all points are read in as the peneplain
pene = ones(size(S.ix));

%begin inversion sequence
 A_built=0
if A_built==0

A_built=0.
%specify expected mean value
u_mean=0.5; %was 100

%specify time step and total time
t_total=2.
t_total=t_total*1000.

lower_bound=0.
upper_bound=2.

%Grid size
%grid_size=1000.

%pi
pi=atan(1.)*4.

disp ('reading in data')

%sort the previously extracted river data ready for the inversion sequence
dx = d;
id = give;
receiver_id = rec;
ps = pene;


lon=x;
lat=y;
number_of_nodes=length(id(:,1));
outlet=id;
ps=double(ps);
n_pene=sum(ps);
receiver=receiver_id;
elev=z0;   


tau_01=chi;
x_grid=zeros([n_pene,1]);
x_grid=zeros([n_pene,1]);
x_grid=zeros([n_pene,1]);
for i = 1:number_of_nodes
 if receiver(i) == 0
        receiver(i) = id(i);
 end
end
k=0;
ij=0;
n_pene=0

for i=1:number_of_nodes
    k=k+1;
    ij=ij+1;
    id2(i)=i;
   
     if id(i)==receiver(i) 
         outlet(i)=1;
         %k=k-1;
         
     end
     
       id_n(k)=i;
       id_n2(k)=i;
       id_row(i)=k;
      %if ps(i)==1;
           n_pene=n_pene+1;
           ipene_row(n_pene)=k;
           ipene(k)=1;
      %end
end
n=k-1;
nnode=ij-1;


disp ('Data read in')
x1=min(x);
x2=max(x);
y1=min(y);
y2=max(y);

x=x-x1;
y=y-y1;
x2=x2-x1;
y2=y2-y1;

nx_grid=1 + floor((x2-0)/grid_size)
ny_grid=1 + floor((y2-0)/grid_size)
nn_grid=nx_grid*ny_grid;

x_grid=zeros([nx_grid,ny_grid]);
y_grid=zeros([nx_grid,ny_grid]);


%BUILD GRID - used for solution
for j=1:ny_grid
    for i=1:nx_grid
       x_grid(i,j)=grid_size*(i-1);
       y_grid(i,j)=grid_size*(j-1);
    end
end

lon_grid=zeros([nx_grid,ny_grid]);
lat_grid=zeros([nx_grid,ny_grid]);

%%PARAMETERS TO CHANGE%%%%%%%%%%%%%%%%%
%specify fluvial erosion rate coefficient
fluv_k=7.;

%Increasing alpha forces smoother spatial solutions
%i.e. block uplift scenarios
%alpha=10.

%Increasing alpha_t forces smoother temporal solutions
%i.e. steady state scenarios
alpha_t=10.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for ifluv_k=1

fluv_k=1

tau=tau_01/fluv_k;
t_total=2*max(tau);
deltat=t_total

m_max=1 
            
%%% BUILD MODEL MATRIX AND DAMPING MATRIX!!!%%%%%%

%First loop through to determine structure

fprintf(1,'Calculating size of A \n');

filled_row=zeros([1+(n*m_max),1]);
filled_row_down=zeros([1+(n*m_max),1]);
filled_row_incision=zeros([1+(n*m_max),1]);
filled_row_grid=zeros([1+(nn_grid*m_max),1]);
count=0;

number_of_entries_in_A=0;
for irow=1:n_pene
    irow2=ipene_row(irow);
    if ipene(irow2)==1 
        filled_row(:)=0.;
        irow1=irow; 
      i=id_n(irow2);
      filled_row=row_fill(irow2,i,m_max,id_row,id_n,receiver,id,tau,deltat,n); 
      filled_row((n*m_max)+1)=1.;
    end
    
    [row,v] = find(filled_row);
    filled_row_grid(:)=0;
    for j=1:n*m_max
        if filled_row(j)>0
            jk = floor((j-1)/(m_max)) +1;
            jk=id_n(jk);
           
            i1=floor(x(jk)/grid_size)+1;
            j1=floor(y(jk)/grid_size)+1;
            spatial_index=i1+(j1-1)*nx_grid;
            temporal_index=mod(j,m_max);
            if temporal_index==0
                temporal_index=m_max;
            end
            full_index=temporal_index+((spatial_index-1)*m_max);
            filled_row_grid(full_index)=filled_row(j);          
        end
    end
    filled_row_grid(1+(nn_grid*m_max))=1.;
    number_of_entries_in_A=number_of_entries_in_A+nnz(filled_row_grid);
end

fprintf(1,'Building A \n');

A_i=zeros([number_of_entries_in_A,1]);
A_j=zeros([number_of_entries_in_A,1]);
A_v=zeros([number_of_entries_in_A,1]);

%Now loop through and fill A matrix

filled_row=zeros([n*m_max,1]);
filled_row_grid=zeros([nn_grid*m_max,1]);

count=0;
number_of_entries_in_A=0;
for irow=1:n_pene
    irow2=ipene_row(irow);
    if ipene(irow2)==1 
        filled_row(:)=0.;
        irow1=irow; 
      i=id_n(irow2);
      filled_row=row_fill(irow2,i,m_max,id_row,id_n,receiver,id,tau,deltat,n); 
      filled_row((n*m_max)+1)=1.;
    end
    
    filled_row_grid(:)=0;
    for j=1:n*m_max
        if filled_row(j)>0
            jk = floor((j-1)/(m_max)) +1;        
            jk=id_n(jk);
           
            i1=floor(x(jk)/grid_size)+1;
            j1=floor(y(jk)/grid_size)+1;
            spatial_index=i1+(j1-1)*nx_grid;
            temporal_index=mod(j,m_max);
            if temporal_index==0
                temporal_index=m_max;
            end
            full_index=temporal_index+((spatial_index-1)*m_max);
            filled_row_grid(full_index)=filled_row_grid(full_index)+filled_row(j);
        end
    end
    filled_row_grid(1+(nn_grid*m_max))=1.;
    for j=1:nn_grid*m_max+1,
        if filled_row_grid(j)>0.
            count=count+1;
            A_i(count)=irow;
            A_j(count)=j;
            A_v(count)=filled_row_grid(j);
        end
    end
end


fprintf(1,'A vectors defined, assembling sparse matrix \n');

A=sparse(A_i,A_j,A_v,n_pene,1+(nn_grid*m_max));


%%NOW Do the same for the full A matrix

number_of_entries_in_A=0;
for irow=1:n
    irow2=irow;
    filled_row(:)=0.;
     irow1=irow; 
      i=id_n(irow2);
      filled_row=row_fill(irow2,i,m_max,id_row,id_n,receiver,id,tau,deltat,n); 
      filled_row((n*m_max)+1)=1.;
    
    [row,v] = find(filled_row);
    filled_row_grid(:)=0;
    for j=1:n*m_max
        if filled_row(j)>0
            jk = floor((j-1)/(m_max)) +1;
            jk=id_n(jk);
           
            i1=floor(x(jk)/grid_size)+1;
            j1=floor(y(jk)/grid_size)+1;
            spatial_index=i1+(j1-1)*nx_grid;
            temporal_index=mod(j,m_max);
            if temporal_index==0
                temporal_index=m_max;
            end
            full_index=temporal_index+((spatial_index-1)*m_max);
            filled_row_grid(full_index)=filled_row(j);          
        end
    end
    filled_row_grid(1+(nn_grid*m_max))=1.;
    number_of_entries_in_A=number_of_entries_in_A+nnz(filled_row_grid);
end

fprintf(1,'Building A \n');

A_full_i=zeros([number_of_entries_in_A,1]);
A_full_j=zeros([number_of_entries_in_A,1]);
A_full_v=zeros([number_of_entries_in_A,1]);

%Now loop through and fill A matrix

filled_row=zeros([n*m_max,1]);
filled_row_grid=zeros([nn_grid*m_max,1]);

count=0;
number_of_entries_in_A=0;
for irow=1:n
    irow2=irow;
    filled_row(:)=0.;
    irow1=irow; 
    i=id_n(irow2);
    filled_row=row_fill(irow2,i,m_max,id_row,id_n,receiver,id,tau,deltat,n); 
    filled_row((n*m_max)+1)=1.;
    
    filled_row_grid(:)=0;
    for j=1:n*m_max
        if filled_row(j)>0
            jk = floor((j-1)/(m_max)) +1;        
            jk=id_n(jk);
           
            i1=floor(x(jk)/grid_size)+1;
            j1=floor(y(jk)/grid_size)+1;
            spatial_index=i1+(j1-1)*nx_grid;
            temporal_index=mod(j,m_max);
            if temporal_index==0
                temporal_index=m_max;
            end
            full_index=temporal_index+((spatial_index-1)*m_max);
            filled_row_grid(full_index)=filled_row_grid(full_index)+filled_row(j);
        end
    end
    filled_row_grid(1+(nn_grid*m_max))=1.;
    for j=1:nn_grid*m_max+1,
        if filled_row_grid(j)>0.
            count=count+1;
            A_full_i(count)=irow;
            A_full_j(count)=j;
            A_full_v(count)=filled_row_grid(j);
        end
    end
end

fprintf(1,'A vectors defined, assembling sparse matrix \n');

A_full=sparse(A_full_i,A_full_j,A_full_v,n,1+(nn_grid*m_max));

fprintf(1,'A matrix built, calculating size of W \n');

end

%for ialpha=r

%ialpha=10
%alpha=1000;  %%%%
    
%Now Build Weighting matrix

%Time could be saved by exploiting symmetry and repetition

count=0
%for time=1:m_max
  for j=1:ny_grid
    for i=1:nx_grid
       weight=4;
       %Look left
       i1=i-1;
       j1=j;
       if i1<1
           weight=weight-1;
       else
           count=count+1;
       end  
       %Look right
       i2=i+1;
       j2=j;
       if i2>nx_grid
           weight=weight-1;
       else
           count=count+1;
       end  
       %Look down
       i3=i;
       j3=j-1;
       if j3<1
           weight=weight-1;
       else
           count=count+1;
       end  
       %Look up
       i4=i;
       j4=j+1;
       if j4>ny_grid
           weight=weight-1;
       else
           count=count+1;
       end  
       count=count+1;
    end
  end
%end

% The additional 1 here is the for the BL parameter
count=count+1

fprintf(1,'Size of W calculated, building W vectors \n');

W_i=zeros([count,1]);
W_j=zeros([count,1]);
W_v=zeros([count,1]);

count=0
for itime=1:m_max
    for j=1:ny_grid
        for i=1:nx_grid
            
            spatial_index=i+(j-1)*nx_grid;
            temporal_index=itime;
            row_index=temporal_index+((spatial_index-1)*m_max);
            
            weight=4;
            %Look left
            i1=i-1;
            j1=j;
            if i1<1
                weight=weight-1;
            else
                count=count+1;
                spatial_index=i1+(j1-1)*nx_grid;
                col_index=temporal_index+((spatial_index-1)*m_max);  
                W_i(count)=row_index;
                W_j(count)=col_index;
                W_v(count)=-1;
            end
            %Look right
            i2=i+1;
            j2=j;
            if i2>nx_grid
                weight=weight-1;
            else
                count=count+1;
                spatial_index=i2+(j2-1)*nx_grid;
                col_index=temporal_index+((spatial_index-1)*m_max);
                W_i(count)=row_index;
                W_j(count)=col_index;
                W_v(count)=-1;
                
            end
            %Look down
            i3=i;
            j3=j-1;
            if j3<1
                weight=weight-1;
            else
                count=count+1;
                spatial_index=i3+(j3-1)*nx_grid;
                col_index=temporal_index+((spatial_index-1)*m_max);
                W_i(count)=row_index;
                W_j(count)=col_index;
                W_v(count)=-1;
            end
            %Look up
            i4=i;
            j4=j+1;
            if j4>ny_grid
                weight=weight-1;
            else
                count=count+1;
                spatial_index=i4+(j4-1)*nx_grid;
                col_index=temporal_index+((spatial_index-1)*m_max);
                W_i(count)=row_index;
                W_j(count)=col_index;
                W_v(count)=-1;
            end
            count=count+1;
            spatial_index=i+(j-1)*nx_grid;
            col_index=temporal_index+((spatial_index-1)*m_max);
            W_i(count)=row_index;
            W_j(count)=col_index;
            W_v(count)=weight;
        end
    end
end

fprintf(1,'Assembling W \n');

%Multiply W_v by alpha where alpha controls the damping
W_v(1:count)=W_v(1:count)*alpha;

count=count+1;
W_i(count)=(nx_grid*ny_grid)+1;
W_j(count)=(nx_grid*ny_grid)+1;
W_v(count)=0.;

%Build Wd
W=sparse(W_i,W_j,W_v,(nx_grid*ny_grid)+1,(nx_grid*ny_grid)+1);

%Build Data weighting matrix

Wd_i=zeros([n,1]);
Wd_j=zeros([n,1]);
Wd_v=zeros([n,1]);

Wd_i=zeros([(n_pene+1)+(1*(nn_grid*m_max)),1]);
Wd_j=zeros([(n_pene+1)+(1*(nn_grid*m_max)),1]);
Wd_v=zeros([(n_pene+1)+(1*(nn_grid*m_max)),1]);

%data weighting
k=1;
for i=1:n_pene
  Wd_i(i)=k;
  Wd_j(i)=k;
  Wd_v(i)=1;
k=k+1;
end

%Idenity matrix for the spatial weighting matrix
for i=1:nn_grid*m_max
Wd_i(k)=k;
Wd_j(k)=k;
Wd_v(k)=1.;
k=k+1;
end    

%BL_weighting again in the spatial weighing matrix
Wd_j(k)=k;
Wd_i(k)=k;
Wd_v(k)=1.;
k=k+1;  

%Build Wd
Wd=sparse(Wd_i,Wd_j,Wd_v);

fprintf(1,'W built, combining matrices \n');

elev_padded=zeros([(n_pene)+1*(1+(nn_grid*m_max)),1]);

for irow=1:n_pene
irow2=ipene_row(irow);
elev_padded(irow)=elev_padded(irow)+elev(irow2);
end

Model_damp=[A;W];
Model_damp=Wd*Model_damp;
%Model_damp=[Wd*A;W;Wt];

elev_padded=Wd*elev_padded;

u_dot=zeros([1+(nn_grid*m_max),1]);
u_pr= zeros([1+(nn_grid*m_max),1]);

u_pr(:)=u_mean;
u_pr(1+(nn_grid*m_max))=0;


%lower and upper bounds on solution
lb=zeros([1+(nn_grid*m_max),1]);
ub=zeros([1+(nn_grid*m_max),1]);

lb(1:nn_grid*m_max)=lower_bound;
ub(1:nn_grid*m_max)=upper_bound;

lb(1+(nn_grid*m_max))=0.000;
ub(1+(nn_grid*m_max))=0.001;

fprintf(1,'Solving system');

u_dot=lsqlin(Model_damp,elev_padded,[],[],[],[],lb,ub,u_pr);

Z=reshape(u_dot(1:nx_grid*ny_grid),nx_grid,ny_grid);
Z = Z*A0^mn

%extract grid coordinates
x_grid_lon = x_grid + x1;
y_grid_lat = y_grid + y1;

%write out the x, y, surface uplift and channel steepness matrices
writematrix(x_grid_lon, 'x_grid.txt')
writematrix(y_grid_lat, 'y_grid.txt')
writematrix(Z, 'ustar.txt')

inversion.xgrid = x_grid_lon;
inversion.ygrid = y_grid_lat;
inversion.ksn = Z;


end