function [rdh_image, image_recovery_check, iteration_max, EC_list, LM_size_list]=mbp(image,iteration_max)
switch nargin
    case 1
        iteration_max = 300;
end

rng(0);
%image = image(5:end-5,5:end-5);
image_size = size(image);

P_s_list=zeros(1,iteration_max);
P_c_list=zeros(1,iteration_max);
EC_list=zeros(1,iteration_max);
LM_size_list=zeros(1,iteration_max);
ref_image_hor = zeros(image_size(1)*image_size(2),iteration_max);
ref_image_hor(:,1) = reshape(image,image_size(1)*image_size(2),1);

image_hor=ref_image_hor(:,1);
original_brightness = mean(ref_image_hor(:,1));
P_s=0;
P_c=0;
iteration=0;
tic
while true
    %Adaptive peak selection
    P_s_previous=P_s;
    P_c_previous=P_c;
    [P_s, P_c]=adaptive_peak_selection(original_brightness,image_hor);
    
    %Direction of histogram shifting
    if P_s < P_c %RHS
        d = 1;
    else %LHS
        d = -1;
    end
    
    %Record location map
    LM=double(image_hor(image_hor==P_c | image_hor==P_c-d)==P_c-d);
    
    %Embedding capacity + stop condition check
    H_P_s=sum(image_hor==P_s);
    
    if H_P_s-sum(image_hor(1:16)==P_s) < length(LM)+32 | iteration == iteration_max %Stop condition reached
        %no need to update iteration
        P_s=P_s_previous;
        P_c=P_c_previous;
        P_s_previous=P_s_list(iteration);
        P_c_previous=P_c_list(iteration);
        first_16_pixels=ref_image_hor(1:16,iteration);
        original_16_lsb=mod(first_16_pixels,2);
        
        if P_s < P_c %RHS
            d = 1;
        else %LHS
            d = -1;
        end
        
        
        %Exclude first 16 pixels from histogram shifting
        image_hor = ref_image_hor(17:end,iteration);
        H_P_s=sum(image_hor==P_s);
        LM=(image_hor(image_hor==P_c | image_hor==P_c-d)==P_c-d);
        payload =randi([0,1],H_P_s-length(LM)-32,1);
        message=[LM ; de2bi(P_s_previous,8)'; de2bi(P_c_previous,8)';original_16_lsb;payload];
        EC_list(iteration)=H_P_s-length(LM)-32;
        LM_size_list(iteration)=length(LM);
        message_whole=zeros(length(image_hor),1);
        message_whole(image_hor==P_s)=message;
        
        %Combine P_c with its neighbor
        image_hor(image_hor==P_c-d)=image_hor(image_hor==P_c-d)+d;
        
        %Shift P_s's neighbors towards P_c
        if d == 1
            image_hor(image_hor > P_s & image_hor < P_c)=image_hor(image_hor > P_s & image_hor < P_c)+d; %RHS
        else
            image_hor(image_hor < P_s & image_hor > P_c)=image_hor(image_hor < P_s & image_hor > P_c)+d; %LHS
        end
        %Embed P_s
        image_hor(image_hor==P_s & message_whole)=image_hor(image_hor==P_s & message_whole)+d;
        
        %Append back the first 16 pixels and replace 16 lsbs with P_s and P_c
        image_hor=[bitxor(bitxor(first_16_pixels,mod(first_16_pixels,2)),[de2bi(P_s,8)'; de2bi(P_c,8)']) ;image_hor];
        iteration_max = iteration;
        rdh_image=reshape(image_hor,image_size(1),image_size(2));
        
        EC_list(iteration_max+1:end)=[];
        LM_size_list(iteration_max+1:end)=[];
        ref_image_hor(:,iteration_max+1:end) = [];
        break
        
    else
        payload =randi([0,1],H_P_s-length(LM)-16,1);
        message=[LM ; de2bi(P_s_previous,8)'; de2bi(P_c_previous,8)';payload];
        iteration=iteration+1;
        EC_list(iteration)=H_P_s-length(LM)-16;
        LM_size_list(iteration)=length(LM);
        
        ref_image_hor(:,iteration)=image_hor;
        P_s_list(iteration)=P_s_previous;
        P_c_list(iteration)=P_c_previous;
    end
    message_whole=zeros(length(image_hor),1);
    message_whole(image_hor==P_s)=message;
    
    
    %Combine P_c with its neighbor
    image_hor(image_hor==P_c-d)=image_hor(image_hor==P_c-d)+d;
    
    %Shift P_s's neighbors towards P_c
    if d == 1
        image_hor(image_hor > P_s & image_hor < P_c)=image_hor(image_hor > P_s & image_hor < P_c)+d; %RHS
    else
        image_hor(image_hor < P_s & image_hor > P_c)=image_hor(image_hor < P_s & image_hor > P_c)+d; %LHS
    end
    %Embed P_s
    image_hor(image_hor==P_s & message_whole)=image_hor(image_hor==P_s & message_whole)+d;
    
end
disp("Encoding time")
toc
% figure(1);imshow(uint8(image));figure(2);imshow(uint8(reshape(ref_image_hor(:,end),image_size(1),image_size(2))))
% figure(3);histogram(image,256);figure(4);histogram(ref_image_hor(:,iteration),256)

tic
% Reverse Operation
first_16_pixels_rec=bitxor(image_hor(1:16),mod(image_hor(1:16),2));
P_s_rec=bi2de(mod(image_hor(1:8)',2));
P_c_rec=bi2de(mod(image_hor(9:16)',2));
% disp("P_s")
% isequal(P_s,P_s_rec)
% disp("P_c")
% isequal(P_c,P_c_rec)
if P_s_rec < P_c_rec %RHS
    d = 1;
else %LHS
    d = -1;
end

%Undo first iteration
image_hor=image_hor(17:end);
%Extract Payload + side information
message_rec=(image_hor(image_hor==P_s_rec |image_hor==P_s_rec+d)==P_s_rec+d);
LM_size=sum(image_hor==P_c_rec);
LM_rec=message_rec(1:LM_size);
% disp("LM")
% isequal(LM,LM_rec)
P_s_p_rec=bi2de(message_rec(1+LM_size:8+LM_size)');
P_c_p_rec=bi2de(message_rec(9+LM_size:16+LM_size)');
first_16_pixels_rec=bitxor(first_16_pixels_rec,message_rec(17+LM_size:32+LM_size));
payload_rec=message_rec(33+LM_size:end);

%Shift back
if d == 1
    image_hor(image_hor > P_s_rec & image_hor < P_c_rec)=image_hor(image_hor > P_s_rec & image_hor < P_c_rec)-d; %RHS
else
    image_hor(image_hor < P_s_rec & image_hor > P_c_rec)=image_hor(image_hor < P_s_rec & image_hor > P_c_rec)-d; %LHS
end
%Undo location map
image_hor(image_hor==P_c_rec)=image_hor(image_hor==P_c_rec)-d*LM_rec;
image_hor=[first_16_pixels_rec; image_hor];

% disp("image_hor")
% isequal(image_hor,ref_image_hor(:,iteration))

iteration=iteration-1;
%Undo rest of the iteration
P_s_rec=P_s_p_rec;
P_c_rec=P_c_p_rec;

% disp("P_s")
% isequal(P_s_rec,P_s_list(iteration+1))
% disp("P_c")
% isequal(P_c_rec,P_c_list(iteration+1))

while (P_s_rec ~= 0 | P_c_rec ~= 0)
    if P_s_rec < P_c_rec %RHS
        d = 1;
    else %LHS
        d = -1;
    end
    
    %Extract Payload + side information
    message_rec=(image_hor(image_hor==P_s_rec |image_hor==P_s_rec+d)==P_s_rec+d);
    LM_size=sum(image_hor==P_c_rec);
    LM_rec=message_rec(1:LM_size);
    P_s_p_rec=bi2de(message_rec(1+LM_size:8+LM_size)');
    P_c_p_rec=bi2de(message_rec(9+LM_size:16+LM_size)');
    payload_rec=message_rec(17+LM_size:end);
    
    %Shift back
    if d == 1
        image_hor(image_hor > P_s_rec & image_hor < P_c_rec)=image_hor(image_hor > P_s_rec & image_hor < P_c_rec)-d; %RHS
    else
        image_hor(image_hor < P_s_rec & image_hor > P_c_rec)=image_hor(image_hor < P_s_rec & image_hor > P_c_rec)-d; %LHS
    end
    %Undo location map
    image_hor(image_hor==P_c_rec)=image_hor(image_hor==P_c_rec)-d*LM_rec;
    P_s_rec=P_s_p_rec;
    P_c_rec=P_c_p_rec;
    
    iteration=iteration-1;
    %     disp("image_hor")
    %     isequal(image_hor,ref_image_hor(:,iteration+1))
    %     [P_s_rec P_c_rec]
end
image_recovery_check=isequal(reshape(image_hor,image_size(1),image_size(2)),image);
disp("Decoding time")
toc
end

function P_c=find_P_c(table,P_s)
combine_table = [table(1:end-1,1)+table(2:end,1) table(1:end-1,2)];
[sort_combine_table, sort_combine_table_index]=sort(combine_table(:,1));
sort_combine_table = [sort_combine_table combine_table(sort_combine_table_index,2)];
list_P_c=sort_combine_table(sort_combine_table(:,1)==sort_combine_table(1,1),:);

%Update list_P_c based on RHS
list_P_c(list_P_c>P_s)=list_P_c(list_P_c>P_s)+1;

%remove P_s-1, P_s, P_s+1 from list_P_c
list_P_c(list_P_c(:,2)==P_s,:)=[];
list_P_c(list_P_c(:,2)==P_s-1,:)=[];
list_P_c(list_P_c(:,2)==P_s+1,:)=[];
[~,index]=min(abs(list_P_c(:,2)-P_s));
P_c=list_P_c(index,2);
end

function P_c=find_P_c_RHS(table,P_s)
combine_table = [table(1:end-1,1)+table(2:end,1) table(2:end,2)];
%Remove entries which cannot be P_c
combine_table(combine_table(:,2)<P_s+2,:)=[];

[sort_combine_table, sort_combine_table_index]=sort(combine_table(:,1));
sort_combine_table = [sort_combine_table combine_table(sort_combine_table_index,2)];
list_P_c=sort_combine_table(sort_combine_table(:,1)==sort_combine_table(1,1),:);

[~,index]=min(abs(list_P_c(:,2)-P_s));
P_c=list_P_c(index,2);
end

function P_c=find_P_c_LHS(table,P_s)
combine_table = [table(1:end-1,1)+table(2:end,1) table(1:end-1,2)];

%Remove entries which cannot be P_c
combine_table(combine_table(:,2)>P_s-2,:)=[];

[sort_combine_table, sort_combine_table_index]=sort(combine_table(:,1));
sort_combine_table = [sort_combine_table combine_table(sort_combine_table_index,2)];
list_P_c=sort_combine_table(sort_combine_table(:,1)==sort_combine_table(1,1),:);
[~,index]=min(abs(list_P_c(:,2)-P_s));
P_c=list_P_c(index,2);
end

function [P_s, P_c]=adaptive_peak_selection(original_brightness,image_hor)
current_brightness = mean(image_hor);
enhancement_mode =0;
if enhancement_mode == 'maximum'
    RHS_limit=-1000;
    LHS_limit=1000;
else
    RHS_limit=0;
    LHS_limit=0;
end

table = [ zeros(256,1) transpose(0:255)];

for i=1:length(image_hor)
    table(image_hor(i)+1)=table(image_hor(i)+1)+1;
end

if current_brightness - original_brightness < RHS_limit %%RHS
    small_table=table(1:254,:);
    [~,index]=max(small_table(:,1));
    P_s=small_table(index,2);
    P_c=find_P_c_RHS(table,P_s);
elseif current_brightness -original_brightness > LHS_limit %%LHS
    small_table=table(3:end,:);
    [~,index]=max(small_table(:,1));
    P_s=small_table(index,2);
    %Find P_c
    P_c=find_P_c_LHS(table,P_s);
else
    [~,index]=max(table(:,1));
    P_s=table(index,2);
    if P_s < 2 %% RHS
        %Find P_c
        P_c=find_P_c_RHS(table,P_s);
        
    elseif P_s > 253 %% LHS
        %Find P_c
        P_c=find_P_c_LHS(table,P_s);
    else %% RHS or LHS
        %Find P_c
        P_c=find_P_c(table,P_s);
    end
    
end
end
