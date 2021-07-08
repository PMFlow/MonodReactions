function [tht]=theta(p)

theta_neg = 1./(3.3333-p);
theta_pos = 0.3;
neg = p < 0; 
tht = theta_neg.*neg + theta_pos.*~neg;

