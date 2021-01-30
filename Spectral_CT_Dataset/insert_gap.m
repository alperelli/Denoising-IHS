
% INSERT_GAP
function data = insert_gap(data,pixel,volD,pix_n,c)
 
if nargin<5
    c = 0.775;
    if nargin<4
        pix_n=2;
        if nargin<3
            volD=0;
        end
    end
end
if(pix_n==2)
    if volD==0
        line = ones(size(data,1),1,size(data,3));
        insert = [line.*(c*data(:,pixel,:)+(1-c)*data(:,pixel+1,:)) line.*((1-c)*data(:,pixel,:)+c*data(:,pixel+1,:))];
        data = [data(:,1:pixel,:) insert data(:,pixel+1:end,:)];
    else
        line = ones(1,size(data,2),size(data,3));
        insert = [line.*(c*data(pixel,:,:)+(1-c)*data(pixel+1,:,:)); line.*((1-c)*data(pixel,:,:)+c*data(pixel+1,:,:))];
        data = [data(1:pixel,:,:); insert; data(pixel+1:end,:,:)];
    end
else if(pix_n==4)
        if volD==0
            line = ones(size(data,1),1,size(data,3));
            insert = [line.*(0.75*data(:,pixel,:)+0.25*data(:,pixel+1,:)) line.*(0.6*data(:,pixel,:)+0.4*data(:,pixel+1,:)) line.*(0.4*data(:,pixel,:)+0.6*data(:,pixel+1,:)) line.*(0.25*data(:,pixel,:)+0.75*data(:,pixel+1,:))];
            data = [data(:,1:pixel,:) insert data(:,pixel+1:end,:)];
        else
            line = ones(1,size(data,2),size(data,3));
            insert = [line.*(0.75*data(:,pixel,:)+0.25*data(:,pixel+1,:)); line.*(0.6*data(:,pixel,:)+0.4*data(:,pixel+1,:)); line.*(0.4*data(:,pixel,:)+0.6*data(:,pixel+1,:)); line.*(0.25*data(:,pixel,:)+0.75*data(:,pixel+1,:))];
            data = [data(1:pixel,:,:); insert; data(pixel+1:end,:,:)];
        end
    end
end
