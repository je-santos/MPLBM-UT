% synthetic geometries

global im_save

im_dir  = 'synthetic_volumes';
im_save = 'binary_volumes_syn';

%% Finneypack operations
num=1;
f1=fopen('synthetic_volumes/spheres_a10_dx0.04_n500_segmented_unsigned_char.raw','r');
d_size = 500;
im=fread(f1,d_size^3,'uint8=>uint8');
im=reshape(im,d_size,d_size,d_size);


SE = strel('sphere',1); %mask
im_dilated = uint8(1-imdilate(~im,SE)); %pore space dilation
save_files(im_dilated, num, 0)

for i=1:8
    save_files(im, num, i)
    im = imdilate(im,SE); %pore space erosion
end

%% fractures
num = 2;
for i=0:4
   im = load(['synthetic_volumes\solid_fracture_6' num2str(i) '.mat']); 
   im = im.domain;
   save_files(im, num, i)
   im = load(['synthetic_volumes\solid_fracture_7' num2str(i) '.mat']); 
   im = im.domain;
   save_files(im, num, i+10)
   im = load(['synthetic_volumes\solid_fracture_8' num2str(i) '.mat']); 
   im = im.domain;
   save_files(im, num, i+20)
   im = load(['synthetic_volumes\solid_fracture_9' num2str(i) '.mat']); 
   im = im.domain;
   save_files(im, num, i+30)
end

%% process-based packs
num = 3;
phis = [10,15,20,25,30,35];
i = 1;
pack_loc = 'synthetic_volumes/SpherePacks_diff_realizations/porosity_0';

for pack_por=phis
    for it=1:3
       im = textread([pack_loc num2str(pack_por) ...
                                    '/Str_' num2str(it) '/S_INPUT.dat']);
       im = reshape(im,256,256,258);
       im(:,:,1)=[];im(:,:,end)=[];
       save_files(im, num, i)
       i = i+1;
    end
end


%% shales
num = 4;
im = load('synthetic_volumes/shale_20_4020');
im = im.shale1;
save_files(im, num, 0);

num = 4;
im = load('synthetic_volumes/shale_30_4030');
im = im.shale2;
save_files(im, num, 1);



%% the weird fractures
num = 5;
im = load('synthetic_volumes/solid_fracture_ANU_19010');
im = im.domain;
save_files(im, num, 0);

num = 5;
im = load('synthetic_volumes/solid_fracturedSP_19011');
im = im.domain;
im(:,1,:  ) = im(:,2,    :  );
im(:,end,:) = im(:,end-1,:  );
im(1,:,:  ) = im(2,:,    :  );
im(end,:,:) = im(end-1,:,:  );
save_files(im, num, 1);

num = 5;
im = load('synthetic_volumes/solid_fracturedSP2_19111');
im = im.domain;
im(:,:,1  ) = im(:,:,2);
im(:,:,end) = im(:,:,end-1);
im(:,1,:  ) = im(:,2,    :  );
im(:,end,:) = im(:,end-1,:  );
save_files(im, num, 2);

num = 5;
im = load('synthetic_volumes/solid_proppedSP_19012');
im = im.domain;
im(1,:,:  ) = im(2,:,    :  );
im(end,:,:) = im(end-1,:,:  );
im(:,1,:  ) = im(:,2,    :  );
im(:,end,:) = im(:,end-1,:  );
save_files(im, num, 3);


%% Yings pack
% small to large
num = 6;
im = load('synthetic_volumes/Vol_3005');
im = im.Vol;
im(:,:,1)=[]; im(:,:,end)=[];
save_files(im, num, 0);

% large to small
num = 6;
im = load('synthetic_volumes/Vol_3006');
im = im.Vol;
im(:,:,1)=[]; im(:,:,end)=[];
save_files(im, num, 1);

% random?
num = 7;
im = load('synthetic_volumes/Vol_3007');
im = im.Vol;
im(:,:,1)=[]; im(:,:,end)=[];
save_files(im, num, 0);

% random?
num = 7;
im = load('synthetic_volumes/Vol_4152');
im = im.Vol;
im(:,:,1)=[]; im(:,:,end)=[];
save_files(im, num, 1);

% random?
num = 7;
im = load('synthetic_volumes/Vol_4201');
im = im.Vol;
im(:,:,1)=[]; im(:,:,end)=[];
save_files(im, num, 2);


% real frac
num = 8;
im = load('synthetic_volumes/Vol_Frac_19006');
im = im.domain;
save_files(im, num, 0);

% fake cracs
num = 9;
for i=1:4
    im = load(['synthetic_volumes/Vol_Frac_300' num2str(i) ]);
    im = im.Vol_Frac;
    im(:,:,1)=[]; im(:,:,end)=[];
    save_files(im, num, i);
end

num = 10;
for i=1:3
    f = fopen(['synthetic_volumes/vuggy' num2str(i) '_256x256x256x8bit.raw'  ]);
    im = reshape( (~fread(f))*1.0, [256,256,256]);
    save_files(im, num, i);
end
