function filled_row = row_fill(irow1,i,m_max,id_row,id_n,receiver,id,tau,deltat,n)

irow=irow1;

q=i;
GlobalTime=tau(q);
itime=1;
LocalTime=deltat;
%irec=receiver(id(i)); 
irec = id(receiver(i));
filled_row=zeros([n*m_max,1]);
%filled_row=zeros([m_max,1]);
while GlobalTime>0
  i_column=((m_max*(irow-1))+(m_max+1-itime));
  row2fill = id_row(receiver(id(q)));
  if row2fill==0 
    NextTau=0.;
  else
    irec=id_n(row2fill);
    NextTau=tau(irec);
  end
  if (GlobalTime-NextTau)>=LocalTime 
    filled_row(i_column)=LocalTime;
    GlobalTime=GlobalTime-LocalTime;
    itime=itime+1;
    LocalTime=deltat;
  else 
    filled_row(i_column)=GlobalTime-NextTau;
    LocalTime=LocalTime-(GlobalTime-NextTau);
    GlobalTime=NextTau;
    q=irec;
    irow=row2fill;
  end
end

end
