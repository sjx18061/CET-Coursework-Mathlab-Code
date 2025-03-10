%================
% constant value:
%================
sx = 0; sy = 4;    		% point s
px = 10; py = 3;   		% point p
qx = 22; qy = 3;   		% point q
p1 = 1;            		% permittivity 1
p2 = 3;            		% permittivity 2
c = 3 * 10.^8;     		% speed of light
j = 1j;            		% imaginary unit
pq_x = px : .0001 : qx; 	% from p to q

% 2.4GHz
f_2_4 = 2.4 * 10^9;            		% 2.4GHz frequency
lambda_2_4 = c ./ f_2_4;        	% calculate wave length when f = 2.4GHz
beta_2_4 = (2 * pi) ./ lambda_2_4;  	% calculate wave number when f = 2.4GHz
pf_2_4 = lambda_2_4 ./ (4 * pi);   	% calculate the phase factor when f = 2.4GHz

% 5.8GHz
f_5_8 = 5.8 * 10^9;            		% 5.8GHz frequency
lambda_5_8 = c ./ f_5_8;        	% calculate wave length when f = 5.8GHz
beta_5_8 = (2 * pi) ./ lambda_5_8;  	% calculate wave number when f = 5.8GHz
pf_5_8 = lambda_5_8 ./ (4 * pi);    	% calculate the phase factor when f = 5.8GHz

% calculate coordinate x of each order
ot = 1 : 1 : 5;					% order for top
ob = 1 : 1 : 5;					% order for bottom
not = 1 : 1 : 3;				% number of terms in order with odd number(top)
net = 1 : 1 : 3;				% number of terms in order with even number(top)
nob = 1 : 1 : 3;				% number of terms in order with odd number(bottom)
neb = 1 : 1 : 3;				% number of terms in order with even number(bottom)
ix = 0;						% constant coordinate x of last image point of each order

iy_t_o = @(ot,not) ((2 .* ot) - (not - 1)) .* sy;			% coordinate y of last image point of odd number of order(top)
iy_t_e = @(ot,net) -1 .* (2 .* (ot - 1) - (net - 1)) .* sy;		% coordinate y of last image point of even number of order(top)
iy_b_o = @(ob,nob) -1 .* ((2 .* ob) - nob) .* sy;			% coordinate y of last image point of odd number of order(bottom)
iy_b_e = @(ob,neb) (2 .* (ob + 1) - neb) .* sy;				% coordinate y of last image point of even number of order(bottom)


%=====================
% calculate direct ray
%=====================
dsp = sqrt(((pq_x - sx) .^2) + ((py - sy) .^2));      	% distance between point s and point p

E_i_2_4 = (1 ./ dsp) .*exp((-j) .* beta_2_4 .* dsp); 	% calculate the electric field of direct ray(2.4GHz)
E_i_5_8 = (1 ./ dsp) .*exp((-j) .* beta_5_8 .* dsp); 	% calculate the electric field of direct ray(5.8GHz)

% change to dB
Ef_g_2_4 = 20 .* log10(abs(E_i_2_4));
Ef_g_5_8 = 20 .* log10(abs(E_i_5_8));                              

%===========================
% calculate first order(top)
%===========================
% coordinate last image point = [0,8]                            
m_ft = (iy_t_o(1,1) - 3) ./ (0 - 10);                    	% gradient of last image point, point last reflection angle and point p
b_ft = py - (m_ft .* pq_x);                  			% calculate c
rx_ft = (6 - b_ft) ./ m_ft;            				% calculate the x of last reflection point [rx,6], ry = 6

drp_ft = sqrt(((0 - pq_x) .^2) + ((iy_t_o(1,1) - py) .^2));    	% calculate total reflection path of first order(top)

d_rs_x_ft = abs(rx_ft - sx);               			% distance x-axis between reflected point and point s
d_ir_y_ft = abs(iy_t_o(1,1) - 6);                        	% distance y-axis between last image point and reflection point
angle_i_ft = atand(d_rs_x_ft ./ d_ir_y_ft);    			% calculate angle of incidence
angle_t_ft = asind((1 ./ sqrt(3)) .* sind(angle_i_ft)); 	% calculate angle of refracted

% calculate the refraction coefficient 
r_c_ft = (cosd(angle_i_ft) - (sqrt(p2 ./ p1) .* cosd(angle_t_ft))) ./ (cosd(angle_i_ft) + (sqrt(p2 ./ p1) .* cosd(angle_t_ft)));   

r_c_t_ft = (r_c_ft) .^1;                               				% calculate total refraction coefficient

E_r_2_4_ft = (1 ./ drp_ft) .*exp((-j) * beta_2_4 .* drp_ft) .* r_c_t_ft;    	% calculate electric field of reflection ray(2.4GHz)
E_r_5_8_ft = (1 ./ drp_ft) .*exp((-j) * beta_5_8 .* drp_ft) .* r_c_t_ft;    	% calculate electric field of reflection ray(5.8GHz)

%==============================
% calculate first order(bottom)
%==============================
% coordinate last image point = [0,-4] 
m_fb = (iy_b_o(1,1) - 3) ./ (0 - 10);                  		% gradient of last image point, point last reflection angle and point p
b_fb = py - (m_fb .* pq_x);                  			% calculate c
rx_fb = (0 - b_fb) ./ m_fb;            				% calculate the x of last reflection point [rx,0], ry = 0

drp_fb = sqrt(((0 - pq_x) .^2) + ((iy_b_o(1,1) - py) .^2));    % calculate total reflection path of first order(bottom)

d_rs_x_fb = abs(rx_fb - sx);               			% distance x-axis between reflected point and point s
d_ir_y_fb = abs(iy_b_o(1,1) - 0);                       	% distance y-axis between last image point and reflection point
angle_i_fb = atand(d_rs_x_fb ./ d_ir_y_fb);    			% calculate angle of incidence
angle_t_fb = asind((1 ./ sqrt(3)) .* sind(angle_i_fb)); 	% calculate angle of refracted

% calculate the refraction coefficient
r_c_fb = (cosd(angle_i_fb) - (sqrt(p2 ./ p1) .* cosd(angle_t_fb))) ./ (cosd(angle_i_fb) + (sqrt(p2 ./ p1) .* cosd(angle_t_fb)));   

r_c_t_fb = (r_c_fb) .^1;                             				% calculate total refraction coefficient

E_r_2_4_fb = (1 ./ drp_fb) .*exp((-j) .* beta_2_4 .* drp_fb) .* r_c_t_fb;    	% calculate electric field of reflection ray(2.4GHz)
E_r_5_8_fb = (1 ./ drp_fb) .*exp((-j) .* beta_5_8 .* drp_fb) .* r_c_t_fb;    	% calculate electric field of reflection ray(5.8GHz)


%==============================================
% calculate total electric field in first order
%==============================================
E_total_2_4_fo = E_i_2_4 + E_r_2_4_ft + E_r_2_4_fb;		% total electric field when f = 2.4GHz
E_total_5_8_fo = E_i_5_8 + E_r_5_8_ft + E_r_5_8_fb;		% total electric field when f = 5.8GHz

% change to dB
Ef_g_2_4_fo = 20 .* log10(abs(E_total_2_4_fo));             		
Ef_g_5_8_fo = 20 .* log10(abs(E_total_5_8_fo));

%===========================
% calculate second order(top)
%===========================
% coordinate last image point = [0,-8] 
m_sot = (iy_t_e(2,1) - 3) ./ (0 - 10);                 		% gradient of last image point, point last reflection angle and point p
b_sot = py - (m_sot .* pq_x);                			% calculate c
rx_sot = (0 - b_sot) ./ m_sot;         				% calculate the x of last reflection point [rx,0], ry = 0

drp_sot = sqrt(((0 - pq_x) .^2) + ((iy_t_e(2,1) - py) .^2));   	% calculate total reflection path

d_rs_x_sot = abs(rx_sot - sx);             			% distance x-axis between reflected point and point s
d_ir_x_sot = abs(iy_t_e(2,1) - 0);                      	% distance y-axis between last image point and reflection point
angle_i_sot = atand(d_rs_x_sot ./ d_ir_x_sot); 			% calculate angle of incidence
angle_t_sot = asind((1 ./ sqrt(3)) .* sind(angle_i_sot)); 	% calculate angle of refracted

% calculate the refraction coefficient
r_c_sot = (cosd(angle_i_sot) - (sqrt(p2 ./ p1) .* cosd(angle_t_sot))) ./ (cosd(angle_i_sot) + (sqrt(p2 ./ p1) .* cosd(angle_t_sot)));   

r_c_t_sot = (r_c_sot) .^2;                            				% calculate total refraction coefficient

E_r_2_4_sot = (1 ./ drp_sot) .*exp((-j) .* beta_2_4 .* drp_sot) .* r_c_t_sot;   % calculate electric field of reflection ray(2.4GHz)
E_r_5_8_sot = (1 ./ drp_sot) .*exp((-j) .* beta_5_8 .* drp_sot) .* r_c_t_sot;   % calculate electric field of reflection ray(5.8GHz)

%==============================
% calculate second order(bottom)
%==============================
% coordinate last image point = [0,16] 
m_sob = (iy_b_e(2,1) - 3) ./ (0 - 10);                 		% gradient of last image point, point last reflection angle and point p
b_sob = py - (m_sob .* pq_x);                			% calculate c
rx_sob = (6 - b_sob) ./ m_sob;         				% calculate the x of last reflection point [rx,6], ry = 6

drp_sob = sqrt(((0 - pq_x) .^2) + ((iy_b_e(2,1) - py) .^2));   	% calculate total reflection path

d_rs_x_sob = abs(rx_sob - sx);             			% distance x-axis between reflected point and point s
d_ir_x_sob = abs(iy_b_e(2,1) - 6);                      	% distance y-axis between last image point and reflection point
angle_i_sob = atand(d_rs_x_sob ./ d_ir_x_sob); 			% calculate angle of incidence
angle_t_sob = asind((1 ./ sqrt(3)) .* sind(angle_i_sob)); 	% calculate angle of refracted

% calculate the refraction coefficient
r_c_sob = (cosd(angle_i_sob) - (sqrt(p2 ./ p1) .* cosd(angle_t_sob))) ./ (cosd(angle_i_sob) + (sqrt(p2 ./ p1) .* cosd(angle_t_sob)));   

r_c_t_sob = (r_c_sob) .^2;                           			% calculate total refraction coefficient

E_r_2_4_sob = (1 ./ drp_sob) .*exp((-j) .* beta_2_4 .* drp_sob) .* r_c_t_sob;   % calculate electric field of reflection ray(2.4GHz)
E_r_5_8_sob = (1 ./ drp_sob) .*exp((-j) .* beta_5_8 .* drp_sob) .* r_c_t_sob;   % calculate electric field of reflection ray(5.8GHz)

%===============================================
% calculate total electric field in second order
%===============================================
E_total_2_4_so = E_total_2_4_fo + E_r_2_4_sot + E_r_2_4_sob;	% total electric field when f = 2.4GHz
E_total_5_8_so = E_total_5_8_fo + E_r_5_8_sot + E_r_5_8_sob;	% total electric field when f = 5.8GHz

% change to dB
Ef_g_2_4_so = 20 * log10(abs(E_total_2_4_so));             		
Ef_g_5_8_so = 20 * log10(abs(E_total_5_8_so));

%===========================
% calculate third order(top)
%===========================
% coordinate last image point = [0,20]
m_tot = (iy_t_o(3,2) - 3) ./ (0 - 10);                 		% gradient of last image point, point last reflection angle and point p
b_tot = py - (m_tot .* pq_x);                			% calculate c
rx_tot = (6 - b_tot) ./ m_tot;         				% calculate the x of last reflection point [rx,6], ry = 6

drp_tot = sqrt(((0 - pq_x) .^2) + ((iy_t_o(3,2) - py) .^2));   % calculate total reflection path

d_rs_x_tot  = abs(rx_tot - sx);             			% distance x-axis between reflected point and point s
d_ir_y_tot  = abs(iy_t_o(3,2) - 6);                      	% distance y-axis between last image point and reflection point
angle_i_tot = atand(d_rs_x_tot ./ d_ir_y_tot); 			% calculate angle of incidence
angle_t_tot = asind((1 ./ sqrt(3)) .* sind(angle_i_tot)); 	% calculate angle of refracted

% calculate the refraction coefficient
r_c_tot = (cosd(angle_i_tot) - (sqrt(p2 ./ p1) .* cosd(angle_t_tot))) ./ (cosd(angle_i_tot) + (sqrt(p2 ./ p1) .* cosd(angle_t_tot)));   

r_c_t_tot = (r_c_tot) .^3;                            				% calculate total refraction coefficient

E_r_2_4_tot = (1 ./ drp_tot) .*exp((-j) .* beta_2_4 * drp_tot) .* r_c_t_tot;    % calculate electric field of reflection ray(2.4GHz)
E_r_5_8_tot = (1 ./ drp_tot) .*exp((-j) .* beta_5_8 * drp_tot) .* r_c_t_tot;    % calculate electric field of reflection ray(5.8GHz)

%==============================
% calculate third order(bottom)
%==============================
% coordinate last image point = [0,-16] 
m_tob = (iy_b_o(3,2) - 3) ./ (0 - 10);                		% gradient of last image point, point last reflection angle and point p
b_tob = py - (m_tob .* pq_x);               	 		% calculate c
rx_tob = (0 - b_tob) ./ m_tob;         				% calculate the x of last reflection point [rx,0], ry = 0

drp_tob = sqrt(((0 - pq_x) .^2) + ((iy_b_o(3,2) - py) .^2));  	% calculate total reflection path

d_rs_x_tob  = abs(rx_tob - sx);             			% distance x-axis between reflected point and point s
d_ir_y_tob  = abs(iy_b_o(3,2) - 0);                     	% distance y-axis between last image point and reflection point
angle_i_tob = atand(d_rs_x_tob ./ d_ir_y_tob); 			% calculate angle of incidence
angle_t_tob = asind((1 ./ sqrt(3)) .* sind(angle_i_tob)); 	% calculate angle of refracted

% calculate the refraction coefficient
r_c_tob = (cosd(angle_i_tob) - (sqrt(p2 ./ p1) .* cosd(angle_t_tob))) ./ (cosd(angle_i_tob) + (sqrt(p2 ./ p1) .* cosd(angle_t_tob)));   

r_c_t_tob = (r_c_tob) .^3;	                            			% calculate total refraction coefficient

E_r_2_4_tob = (1 ./ drp_tob) .*exp((-j) .* beta_2_4 .* drp_tob) .* r_c_t_tob;   % calculate electric field of reflection ray(2.4GHz)
E_r_5_8_tob = (1 ./ drp_tob) .*exp((-j) .* beta_5_8 .* drp_tob) .* r_c_t_tob; 	% calculate electric field of reflection ray(5.8GHz)

%==============================================
% calculate total electric field in third order
%==============================================
E_total_2_4_to = E_total_2_4_so + E_r_2_4_tot + E_r_2_4_tob;	% total electric field when f = 2.4GHz
E_total_5_8_to = E_total_5_8_so + E_r_5_8_tot + E_r_5_8_tob;	% total electric field when f = 5.8GHz

% change to dB
Ef_g_2_4_to = 20 .* log10(abs(E_total_2_4_to));             		
Ef_g_5_8_to = 20 .* log10(abs(E_total_5_8_to)); 

%============================
% calculate fourth order(top)
%============================
% coordinate last image point = [0,-20]; 
m_fot = (iy_t_e(4,2) - 3) ./ (0 - 10);                		% gradient of last image point, point last reflection angle and point p
b_fot = py - (m_fot .* pq_x);                			% calculate c
rx_fot = (0 - b_fot) ./ m_fot;         				% calculate the x of last reflection point [rx,0], ry = 0

drp_fot = sqrt(((0 - pq_x) .^2) + ((iy_t_e(4,2) - py) .^2));   	% calculate total reflection path

d_rs_x_fot = abs(rx_fot - sx);             			% distance x-axis between reflected point and point s
d_ir_y_fot = abs(iy_t_e(4,2) - 0);                     		% distance y-axis between last image point and reflection point
angle_i_fot = atand(d_rs_x_fot ./ d_ir_y_fot); 			% calculate angle of incidence
angle_t_fot = asind((1 ./ sqrt(3)) .* sind(angle_i_fot)); 	% calculate angle of refracted

% calculate the refraction coefficient
r_c_fot = (cosd(angle_i_fot) - (sqrt(p2 ./ p1) .* cosd(angle_t_fot))) ./ (cosd(angle_i_fot) + (sqrt(p2 ./ p1) .* cosd(angle_t_fot)));

r_c_t_fot = (r_c_fot) .^4;	                            			% calculate total refraction coefficient

E_r_2_4_fot = (1 ./ drp_fot) .*exp((-j) .* beta_2_4 .* drp_fot) .* r_c_t_fot;   % calculate electric field of reflection ray(2.4GHz)
E_r_5_8_fot = (1 ./ drp_fot) .*exp((-j) .* beta_5_8 .* drp_fot) .* r_c_t_fot; 	% calculate electric field of reflection ray(5.8GHz)

%===============================
% calculate fourth order(bottom)
%===============================
% coordinate last image point = [0,28]; 
m_fob = (iy_b_e(4,2) - 3) ./ (0 - 10);                 		% gradient of last image point, point last reflection angle and point p
b_fob = py - (m_fob .* pq_x);                			% calculate c
rx_fob = (6 - b_fob) ./ m_fob;         				% calculate the x of last reflection point [rx,6], ry = 6

drp_fob = sqrt(((0 - pq_x) .^2) + ((iy_b_e(4,2) - py) .^2));   	% calculate total reflection path

d_rs_x_fob = abs(rx_fob - sx);            	 		% distance x-axis between reflected point and point s
d_ir_y_fob = abs(iy_b_e(4,2) - 6);                      	% distance y-axis between last image point and reflection point
angle_i_fob = atand(d_rs_x_fob ./ d_ir_y_fob); 			% calculate angle of incidence
angle_t_fob = asind((1 ./ sqrt(3)) .* sind(angle_i_fob)); 	% calculate angle of refracted

% calculate the refraction coefficient
r_c_fob = (cosd(angle_i_fob) - (sqrt(p2 ./ p1) .* cosd(angle_t_fob))) ./ (cosd(angle_i_fob) + (sqrt(p2 ./ p1) .* cosd(angle_t_fob)));   

r_c_t_fob = (r_c_fob) .^4;	                            			% calculate total refraction coefficient

E_r_2_4_fob = (1 ./ drp_fob) .*exp((-j) .* beta_2_4 .* drp_fob) .* r_c_t_fob;   % calculate electric field of reflection ray(2.4GHz)
E_r_5_8_fob = (1 ./ drp_fob) .*exp((-j) .* beta_5_8 .* drp_fob) .* r_c_t_fob;   % calculate electric field of reflection ray(5.8GHz)

%===============================================
% calculate total electric field in fourth order
%===============================================
E_total_2_4_foo = E_total_2_4_to + E_r_2_4_fot + E_r_2_4_fob;	% total electric field when f = 2.4GHz
E_total_5_8_foo = E_total_5_8_to + E_r_5_8_fot + E_r_5_8_fob;	% total electric field when f = 5.8GHz

% change to dB
Ef_g_2_4_foo = 20 .* log10(abs(E_total_2_4_foo));           		
Ef_g_5_8_foo = 20 .* log10(abs(E_total_5_8_foo)); 

%===========================
% calculate fifth order(top)
%===========================
% coordinate last image point = [0,32]; 
m_Fot = (iy_t_o(5,3) - 3) ./ (0 - 10);                 		% gradient of last image point, point last reflection angle and point p
b_Fot = py - (m_Fot .* pq_x);                			% calculate c
rx_Fot = (6 - b_Fot) ./ m_Fot;         				% calculate the x of last reflection point [rx,6], ry = 6

drp_Fot = sqrt(((0 - pq_x) .^2) + ((iy_t_o(5,3) - py) .^2));   % calculate total reflection path

d_rs_x_Fot = abs(rx_Fot - sx);             			% distance x-axis between reflected point and point s
d_ir_y_Fot = abs(iy_t_o(5,3) - 6);                      	% distance y-axis between last image point and reflection point
angle_i_Fot = atand(d_rs_x_Fot ./ d_ir_y_Fot); 			% calculate angle of incidence
angle_t_Fot = asind((1 ./ sqrt(3)) .* sind(angle_i_Fot)); 	% calculate angle of refracted

% calculate the refraction coefficient
r_c_Fot = (cosd(angle_i_Fot) - (sqrt(p2 ./ p1) .* cosd(angle_t_Fot))) ./ (cosd(angle_i_Fot) + (sqrt(p2 ./ p1) .* cosd(angle_t_Fot)));

r_c_t_Fot = (r_c_Fot) .^5;	                            				% calculate total refraction coefficient

E_r_2_4_Fot = (1 ./ drp_Fot) .*exp((-j) .* beta_2_4 .* drp_Fot) .* r_c_t_Fot;   	% calculate electric field of reflection ray(2.4GHz)
E_r_5_8_Fot = (1 ./ drp_Fot) .*exp((-j) .* beta_5_8 .* drp_Fot) .* r_c_t_Fot;   	% calculate electric field of reflection ray(5.8GHz)

%==============================
% calculate fifth order(bottom)
%==============================
% coordinate last image point = [0,-28]; 
m_Fob = (iy_b_o(5,3) - 3) ./ (0 - 10);                		% gradient of last image point, point last reflection angle and point p
b_Fob = py - (m_Fob .* pq_x);                			% calculate c
rx_Fob = (0 - b_Fob) ./ m_Fob;        	 			% calculate the x of last reflection point [rx,0], ry = 0

drp_Fob = sqrt(((0 - pq_x) .^2) + ((iy_b_o(5,3) - py) .^2));   	% calculate total reflection path

d_rs_x_Fob = abs(rx_Fob - sx);             			% distance x-axis between reflected point and point s
d_ir_y_Fob = abs(iy_b_o(5,3) - 0);                     		% distance y-axis between last image point and reflection point
angle_i_Fob = atand(d_rs_x_Fob ./ d_ir_y_Fob); 			% calculate angle of incidence
angle_t_Fob = asind((1 ./ sqrt(3)) .* sind(angle_i_Fob)); 	% calculate angle of refracted

% calculate the refraction coefficient
r_c_Fob = (cosd(angle_i_Fob) - (sqrt(p2 ./ p1) .* cosd(angle_t_Fob))) ./ (cosd(angle_i_Fob) + (sqrt(p2 ./ p1) .* cosd(angle_t_Fob)));

r_c_t_Fob = (r_c_Fob) .^5;	                            			% calculate total refraction coefficient

E_r_2_4_Fob = (1 ./ drp_Fob) .*exp((-j) .* beta_2_4 .* drp_Fob) .* r_c_t_Fob;   % calculate electric field of reflection ray(2.4GHz)
E_r_5_8_Fob = (1 ./ drp_Fob) .*exp((-j) .* beta_5_8 .* drp_Fob) .* r_c_t_Fob;   % calculate electric field of reflection ray(5.8GHz)

%==============================================
% calculate total electric field in fifth order
%==============================================
E_total_2_4_Fo = E_total_2_4_foo + E_r_2_4_Fot + E_r_2_4_Fob;	% total electric field when f = 2.4GHz
E_total_5_8_Fo = E_total_5_8_foo + E_r_5_8_Fot + E_r_5_8_Fob;	% total electric field when f = 5.8GHz

% change to dB
Ef_g_2_4_Fo = 20 .* log10(abs(E_total_2_4_Fo));             		
Ef_g_5_8_Fo = 20 .* log10(abs(E_total_5_8_Fo));

%============
% plot graph
%============
% graph for 2.4GHz
figure(1);
hold on;
plot(pq_x, Ef_g_2_4, 'k-', 'LineWidth', 0.5);    	% direct ray
plot(pq_x, Ef_g_2_4_fo, 'c-', 'LineWidth', 0.5); 	% first order
plot(pq_x, Ef_g_2_4_so, 'g-', 'LineWidth', 0.5);  	% second order
plot(pq_x, Ef_g_2_4_to, 'b-', 'LineWidth', 0.5);  	% third order
plot(pq_x, Ef_g_2_4_foo,'r-', 'LineWidth', 0.5);  	% fourth order
plot(pq_x, Ef_g_2_4_Fo, 'm-', 'LineWidth', 0.5);  	% fifth order
hold off;
xlabel('Line Segment "pq"/m'), ylabel('E-field/dB');
title('Total Electric field in each order when f = 2.4GHz');
legend('Direct Ray', 'First Order', 'Second Order', 'Third Order', 'Fourth Order', 'Fifth Order');
grid on;

% graph for 5.8GHz
figure(2);
hold on;
plot(pq_x, Ef_g_5_8, 'k-', 'LineWidth', 0.5);      	% direct ray
plot(pq_x, Ef_g_5_8_fo, 'c-',  'LineWidth', 0.5);   	% first order
plot(pq_x, Ef_g_5_8_so, 'g-',  'LineWidth', 0.5);   	% second order
plot(pq_x, Ef_g_5_8_to, 'b-',  'LineWidth', 0.5);   	% third order
plot(pq_x, Ef_g_5_8_foo,'r-',  'LineWidth', 0.5);   	% fourth order
plot(pq_x, Ef_g_5_8_Fo, 'm-',  'LineWidth', 0.5);   	% fifth order
hold off;
xlabel('Line Segment "pq"/m'), ylabel('E-field/dB');
title('Total Electric field in each order when f = 5.8GHz');
legend('Direct Ray', 'First Order', 'Second Order', 'Third Order', 'Fourth Order', 'Fifth Order');
grid on;

% graph of the relationship between incidence ray and reflection coeeficient
o = 1 : 1 :10; 						% order 1-10 for plot graph
angle_x = 0 : 1 : 90;					% incident angle 0-90 degree
angle_x_t = asind((1 ./ sqrt(3)) .* sind(angle_x));	% refraction angle

% reflection coeeficient
r_c_x = (cosd(angle_x) - (sqrt(p2 ./ p1) .* cosd(angle_x_t))) ./ (cosd(angle_x) + (sqrt(p2 ./ p1) .* cosd(angle_x_t)));
	
r_c_x_t = @(o) r_c_x .^o;				% total reflection coefficient in each order

figure(3);
hold on;

% 1st order
plot(angle_x, r_c_x_t(1), 'Color', [1, 0.84, 0], 'LineWidth', 0.5, 'DisplayName', '1st Order'); 	
% 2nd order	
plot(angle_x, r_c_x_t(2), 'Color', [1, 0.5, 0.31], 'LineWidth', 0.5, 'DisplayName', '2nd Order');
% 3rd order 	
plot(angle_x, r_c_x_t(3), 'Color', [0.2, 0.8, 0.2], 'LineWidth', 0.5, 'DisplayName', '3rd Order'); 
% 4th order 		
plot(angle_x, r_c_x_t(4), 'Color', [0.5, 0.5, 0], 'LineWidth', 0.5, 'DisplayName', '4th Order');
% 5th order  		
plot(angle_x, r_c_x_t(5), 'Color', [0, 0.5, 0.5], 'LineWidth', 0.5, 'DisplayName', '5th Order');
% 6th order  		
plot(angle_x, r_c_x_t(6), 'Color', [1, 0.55, 0], 'LineWidth', 0.5, 'DisplayName', '6th Order');
% 7th order
plot(angle_x, r_c_x_t(7), 'Color', [1, 0, 1], 'LineWidth', 0.5, 'DisplayName', '7th Order'); 
% 8th order 		
plot(angle_x, r_c_x_t(8), 'Color', [0.98, 0.5, 0.45], 'LineWidth', 0.5, 'DisplayName', '8th Order');
% 9th order
plot(angle_x, r_c_x_t(9), 'Color', [1, 0.75, 0], 'LineWidth', 0.5, 'DisplayName', '9th Order');
% 10th order
plot(angle_x, r_c_x_t(10), 'Color', [0.8, 0.25, 0.33], 'LineWidth', 0.5, 'DisplayName', '10th Order');

hold off;
xlabel('Angle Incident/degree'), ylabel('Total Reflection Coeeficient');
title ('Angle Incident VS Total Reflection Coefficient of Each Order');
legend('show', 'Location', 'best');
grid on;

