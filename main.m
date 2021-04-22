function main
%Image
image=double(imread('Kodak images/Original/kodim01_org.png'));

%Payload
rng(1) %Presets randomness to ensure the results are reproducible
payload_length=15000; %number of bits to be embedded
payload=randi([0,1],payload_length,1);

%% Embedding for maximum contrast
iteration_max=1000; % default maximum iteration_max is 1000, you can increase it as much as you want, but it may take longer time
max_contrast_bypass_mode=0; % 0 = embedes addtional synthetic bits after the specified payload has been embedded to maximize the contrast
[rdh_image, ~, ~, ~,embedding_capacity_left]=mbp(image,payload,iteration_max,max_contrast_bypass_mode);
if embedding_capacity_left < 0
    disp('Failed embedding, try increasing iteration_max') 
else
    disp(['Can embed ' num2str(embedding_capacity_left) ' bits more (estimated)'])
end

%% Recovery check for maximum contrast case
[payload_rec, re_image] = mbp_recovery(rdh_image);
if isequal(re_image,image)
    disp('Original image recovered')
else
    disp('Failed to recover the original image')
end

if isequal(payload_rec,payload)
    disp('Payload recovered')
else
    disp('Failed to recover the payload')
end


%% Embedding only the payload => does not maximize the contrast
iteration_max=1000; % default maximum iteration_max is 1000, you can increase it as much as you want, but it may take longer time
max_contrast_bypass_mode=1; %1=terminates after specified payload has been embedded, does not achieve maximum contrast
[rdh_image_non_max_contrast, ~, ~, ~,embedding_capacity_left_non_max_contrast]=mbp(image,payload,iteration_max,max_contrast_bypass_mode);

if embedding_capacity_left_non_max_contrast ~= -1
    disp('Failed embedding, try increasing iteration_max') 
else
    disp(['Successfully embedded ' num2str(payload_length) ' bits. Cannot determine how many more bits can be embedded since max_contrast_bypass_mode was 1. To determine maximum number of bits embeddable (estimated), run with max_contrast_bypass_mode =0'])
end


%% Recovery check for when only the specified payload has been embedded => does not maximize the contrast
[payload_rec_non_max_contrast, re_image_non_max_contrast] = mbp_recovery(rdh_image_non_max_contrast);
if isequal(re_image_non_max_contrast,image)
    disp('Original image recovered')
else
    disp('Failed to recover the original image')
end

if isequal(payload_rec_non_max_contrast,payload)
    disp('Payload recovered')
else
    disp('Failed to recover the payload')
end


%show the image
close all
figure(1)
imshow(uint8(image))
figure(2)
imshow(uint8(rdh_image))
figure(3)
imshow(uint8(rdh_image_non_max_contrast))














end