# LDC_MGM-matlab
local density clustering matlab code


```
clear
close all
tic
addpath E:\local_density_clustering\model
% 针对数据做局部密度聚类
% 参数设置
para.rhomin=1.5;
para.deltamin=4;
para.v_min=64;
para.rms=0.7;
para.sigma=0.6;
para.gradtmin=0.1;
is_plot = 0;
get_out_mask=1;

work_path = 'high_density';
file_list = dir([work_path '\txt_pix\*.txt']);
file_num = length(file_list);

for i=1:file_num
    txt_name = file_list(i).name;
    out_name = [work_path '\fits\out\s_out_' txt_name(10:12) '.fits'];
    model_name = [work_path '\fits\model\s_model_' txt_name(10:12) '.fits'];
    mask_name = [work_path '\result\mask\mask_' txt_name(10:12) '.fits'];
    Outcat_name = [work_path '\result\txt_mat\Outcat_', txt_name(10:12)];
    
    data = fitsread(out_name);
    outcat_name = [work_path '\result\txt_mat1\outcat_', txt_name(10:12)];
    
    [outcat, ~, mask] = localDenClust2_0_1(data,para,is_plot,get_out_mask);
    
    save(outcat_name,'outcat')    
    save_outcat(Outcat,Outcat_name)
    save_mask(mask,mask_name)
  
end
```
