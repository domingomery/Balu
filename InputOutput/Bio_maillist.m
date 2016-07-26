% Bio_maillist(mymail,mypassword,mails,subject,heads,body,signature)
%
% Toolbox: Balu
%    Send e-mail list
%    
%    Send an e-mail form mymail to mails, with subject, message, head(s) and
%    signature. Ir requires the password of mymail.
%
% mails     = {'peter@gmail.com','rosa@hotmail.com','tomas@yahoo.es'};
% heads     = {'Hi Peter:','Dear Rosa:','Hello Tomas:'};
% message   = {'Please do not forget the next meeting on Monday.'}; 
% subject   = {'Next meeting'}; 
% signature = {'Regards',' ','John',' ','----------------------------------','Prof. John Schmidt','Departmento of Computer Science','University of ...','http://www.usw.edu','----------------------------------'};
% Bio_maillist('johns@gmail.com','johnssszzz12',mails,subject,heads,message,signature);
%
% See also Bio_sendmail.
%
% (c) GRIMA-DCCUC, 2012
% http://grima.ing.puc.cl

function Bio_maillist(mymail,mypassword,mails,subject,heads,body,signature)

n = length(mails);
n1 = length(heads);
if n~=n1
    error('number of heads must be equal to number of mails...');
end
nb = length(body);
if and(nb~=1,nb~=n)
    error('number of bodies must be equal to 1 or equal to number of mails...');
end

ns = length(subject);
if and(ns~=1,ns~=n)
    error('number of subjects must be equal to 1 or equal to number of mails...');
end

if nb==1
    msg = body;
end

if ns==1
    sbj = char(subject);
end

for i=1:n
    if nb>1
        msg = body{i};
    end
    if ns>1
        sbj = char(subject{i});
    end
    Bio_sendmail(mymail,mypassword,mails{i},sbj,[heads{i} ' ' msg ' ' signature]);
end
  
  