function B3f=splineB3(x2,r)
x2=x2-(x2(end)+x2(1))/2;
stepx=x2(2)-x2(1);
ri=ceil(2*r); 

r=r*stepx;
x2c=x2((ceil((end+1)/2)));
x=x2(ceil((end+1)/2)-ri:ceil((end+1)/2)+ri);
d=abs(x-x2c)/r;
B3=zeros(size(x));

    for ix=-ri:ri
           id=ix+ri+1;
        if(d(id)<1)%use the first polynomial  
            B3(id)=(3*d(id).^3-6*d(id).^2+4)/6;  
        else if(d(id)<2)
            B3(id)=(-d(id).^3+6*d(id).^2-12*d(id)+8)/6;
            end
        end
    end
B3f=zeros(size(x2));B3f(ceil((end+1)/2)-ri:ceil((end+1)/2)+ri)=B3;

end

