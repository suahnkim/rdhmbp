function main
%Image
image=double(rgb2gray(imread('Kodak/kodim01.png')));

%Payload
payload_length=25000; %number of bits to be embedded
payload=randi([0,1],payload_length,1);

%Embedding
[rdh_image, ~, ~, ~,embedding_capacity_left]=mbp(image,payload);
if embedding_capacity_left < 0
    disp('Failed embedding')
else
    disp(['Can embed ' num2str(embedding_capacity_left) ' bits more (estimated)'])
end

%Recovery check
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

%show the image
close all
figure(1)
imshow(uint8(image))
figure(2)
imshow(uint8(rdh_image))
end