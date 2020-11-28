function repeatability_demo

fprintf(1,'It may take a while, i.e. 30min \n');
mex c_eoverlap.cxx;

% in windows we need to execute the .ln files with wsl
use_windows = 'wsl -- ';

% exchange for other data sets
prefix_rel = "../data/graf/";
prefix_res = sprintf('%sresults/',prefix_rel);

det_suffix=['siftaf';'haraff';'hesaff';'mseraf';'ibraff';'ebraff';'t-01-1';'t-25-2';'t-50-2';'t-20-2'];
%det_suffix=['t-01-1';'t-25-2';'t-50-2';'t-75-2';'t-20-2'];
det_nb=3;


for i=1:6
image = sprintf('%simg%d',prefix_rel,i);
res = sprintf('%simg%d', prefix_res,i);
detectFeatures(sprintf('%s./h_affine.ln -haraff -sift -i %s.ppm -o %s.%s -thres 1000',use_windows, image,res,det_suffix(2,:)));
detectFeatures(sprintf('%s./h_affine.ln -hesaff -sift -i %s.ppm  -o %s.%s -thres 500',use_windows,image,res,det_suffix(3,:)));
detectFeatures(sprintf('%s./mser.ln  -t 2 -es 2 -i %s.ppm -o %s.%s',use_windows,image,res,det_suffix(4,:)));
detectFeatures(sprintf('%s./ibr.ln %s.ppm %s.%s -scalefactor 1.0',use_windows,image,res,det_suffix(5,:)));
detectFeatures(sprintf('%s./ebr.ln %s.ppm %s.%s',use_windows,image,res,det_suffix(6,:)));
end

figure(1);clf;
grid on;
ylabel('repeatebility %')
xlabel('viewpoint angle');
hold on;
figure(2);clf;
grid on;
ylabel('nb of correspondences')
xlabel('viewpoint angle');
hold on;

mark=['-k>';'-ks';'-k+';'-k*';'-kd';'-kv';'-rx';'-gx';'-bx';'-yx'];
%mark=['-rx';'-gx';'-bx';'-yx';'-m*'];

for d=1:10
seqrepeat=[];
seqcorresp=[];
for i=2:6
  file1=sprintf('%simg1.%s',prefix_res,det_suffix(d,:));
file2=sprintf('%simg%d.%s',prefix_res,i,det_suffix(d,:));
Hom=sprintf('%sH1to%dp',prefix_rel,i);
imf1='img1.ppm';
imf2=sprintf('img%d.ppm',i);
[erro,repeat,corresp, match_score,matches, twi]=repeatability(file1,file2,Hom,imf1,imf2, 1);
seqrepeat=[seqrepeat repeat(4)];
seqcorresp=[seqcorresp corresp(4)];
end
figure(1);  plot([20 30 40 50 60],seqrepeat,mark(d,:));
figure(2);  plot([20 30 40 50 60],seqcorresp,mark(d,:));
end

figure(1);legend(det_suffix(1,:),det_suffix(2,:),det_suffix(3,:),det_suffix(4,:),det_suffix(5,:),det_suffix(6,:),det_suffix(7,:),det_suffix(8,:),det_suffix(9,:),det_suffix(10,:));
axis([10 70 0 100]);
figure(2);legend(det_suffix(1,:),det_suffix(2,:),det_suffix(3,:),det_suffix(4,:),det_suffix(5,:),det_suffix(6,:),det_suffix(7,:),det_suffix(8,:),det_suffix(9,:),det_suffix(10,:));

function detectFeatures(command)
fprintf(1,'Detecting features: %s\n',command);
[status,result] = system(command);
