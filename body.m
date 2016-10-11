classdef body < handle
    % Filename: body.m
    % Author:   Samuel Acuña
    % Date:     11 Oct 2016
    % About:
    % a generic spatial rigid body, including all the attributes of the
    % body, as well as local points on the body.This also
    % estabilishes the local reference frame.
   properties
      ID;   % body ID number
      r; % = [x;y;z] location of body RF from GLOBAL RF
      p; % = [e0;e1;e2;e3] euler parameters of BODY RF
      m; % mass of the part
      J; % inertia tensor of the part
   end
   
   methods
       function body = body(ID,r,p) %constructor function
            body.ID = ID;
            body.r = r;
            body.p = p;    
       end


   end
end 