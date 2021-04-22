function [rdh_image, iteration_max, EC_list, LM_size_list, embedding_capacy_left]=mbp(image, actual_payload,iteration_max, max_contrast_bypass_mode)
if isempty(iteration_max)
        iteration_max = 1000;
end
if isempty(max_contrast_bypass_mode)
    max_contrast_bypass_mode = 0;
end

%Preprocess Payload (length appended)
image_size=size(image);
payload_length_max=2*ceil(log2(image_size(1)*image_size(2)+1));
actual_payload=[de2bi(length(actual_payload),payload_length_max)'; actual_payload];

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
payload_total=[];
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
    
    if H_P_s-sum(image_hor(1:16)==P_s) < length(LM)+32 || iteration == iteration_max || (max_contrast_bypass_mode==1 && length(payload_total) > length(actual_payload)) %Stop condition reached
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
        
        payload_total(end-payload_length_last+1:end)=[];
        if length(payload_total) < length(actual_payload)
            payload_left_over=length(actual_payload)-length(payload_total);
            if payload_left_over < H_P_s-length(LM)-32
                synthetic_payload =randi([0,1],H_P_s-length(LM)-32-payload_left_over,1);
                payload = [actual_payload(length(payload_total)+1:end); synthetic_payload];
            else
                payload = actual_payload(length(payload_total)+1:length(payload_total)+H_P_s-length(LM)-32);
            end
        else
            payload =randi([0,1],H_P_s-length(LM)-32,1);
        end
        
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
        payload_total=[payload_total; payload];
        embedding_capacy_left=length(payload_total)-length(actual_payload);
        if max_contrast_bypass_mode ==1
            embedding_capacy_left=-1;
        end
        break
        
    else
        if length(payload_total) < length(actual_payload)
            payload_left_over=length(actual_payload)-length(payload_total);
            if payload_left_over < H_P_s-length(LM)-16
                synthetic_payload =randi([0,1],H_P_s-length(LM)-16-payload_left_over,1);
                payload = [actual_payload(length(payload_total)+1:end); synthetic_payload];
            else
                payload = actual_payload(length(payload_total)+1:length(payload_total)+H_P_s-length(LM)-16);
            end
        else
            payload =randi([0,1],H_P_s-length(LM)-16,1);
        end
        
        payload_length_last=length(payload);
        payload_total=[payload_total; payload];
            
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
rdh_image=reshape(image_hor,image_size(1),image_size(2));

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
