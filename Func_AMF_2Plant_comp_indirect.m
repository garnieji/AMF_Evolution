
function dy = Func_AMF_2Plant_comp_indirect(t,y)


%3 piante, 
global rH qhp qcm qcp qhm mup mum d alpha...
    beta1 beta2 ap12 ap21 

p1 = y(1);
p2 = y(2);
f1 = y(3);
f2 = y(4);
dy = zeros(4,1);

dy(1) = qhp*rH*p1 + qhp*alpha*f1*(ap21)/(ap21+f2)*p1/(d+p1) - ...
    qcp*beta1*f1*p1 - mup*p1^2;
dy(2) = qhp*rH*p2 + qhp*alpha*f2*(ap12)/(ap12+f1)*p2/(d+p2) - ...
    qcp*beta2*f2*p2 - mup*p2^2;
dy(3) = qcm*beta1*f1*p1 - qhm*alpha*f1*(ap21)/(ap21+f2)*p1/(d+p1) - mum*f1^2;
dy(4) = qcm*beta2*f2*p2 - qhm*alpha*f2*(ap12)/(ap12+f1)*p2/(d+p2) - mum*f2^2;


end