function [payload_rec, re_image]=mbp_recovery(rdh_image)
tic
image_size = size(rdh_image);
image_hor = reshape(rdh_image,image_size(1)*image_size(2),1);
% ref_image_hor(:,1) = reshape(rdh_image,image_size(1)*image_size(2),1);

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

%Undo rest of the iteration
P_s_rec=P_s_p_rec;
P_c_rec=P_c_p_rec;

% disp("P_s")
% isequal(P_s_rec,P_s_list(iteration+1))
% disp("P_c")
% isequal(P_c_rec,P_c_list(iteration+1))

while (P_s_rec ~= 0 || P_c_rec ~= 0)
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
    payload_rec=[ message_rec(17+LM_size:end); payload_rec];
    
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
        %     disp("image_hor")
    %     isequal(image_hor,ref_image_hor(:,iteration+1))
    %     [P_s_rec P_c_rec]
end
re_image=reshape(image_hor,image_size(1),image_size(2));
payload_length_max=2*ceil(log2(image_size(1)*image_size(2)+1)); 
payload_length=bi2de(payload_rec(1:payload_length_max)');
payload_rec(1:payload_length_max)=[];
payload_rec(payload_length+1:end)=[];
toc
end

