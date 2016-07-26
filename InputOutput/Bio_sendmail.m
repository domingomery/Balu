% Bio_sendmail(mymail,mypassword,mailto,subject)
% Bio_sendmail(mymail,mypassword,mailto,subject,message)
% Bio_sendmail(mymail,mypassword,mailto,subject,message,attachment)
%
% Toolbox: Balu
%    Send e-mail
%    
%    Send an e-mail form mymail to mailto, with subject, message (optionsl)
%    and attachment (optional). Ir requires the password of mymail.
%
% Example 1:
%
% x=10;y=20;
% msg1 = sprintf('x=%d',x)
% msg2 = sprintf('y=%d',y)
% Bio_sendmail('johns@gmail.com','johnssszzz12','marys@gmail.com','New results',{'Hi Mary,','here are the results:',msg1,msg2,'John'});
%
%
% Example 2:
%
% Bio_sendmail('johns@gmail.com','johnssszzz12','marys@gmail.com','Output Image','Hi Mary, I am sending you the image. John','image.jpg');
%
% See also sendmail.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl
%

function Bio_sendmail(mymail,mypassword,mailto,subject,message,attachment)

fprintf('Sending e-mail to %s...\n',mailto);
setpref('Internet','E_mail',mymail);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mymail);
setpref('Internet','SMTP_Password',mypassword);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class','javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

if exist('message','var')
    if exist('attachment','var')
        sendmail(mailto,subject,message,attachment);
        fprintf('>> attaching file %s...\n',attachment);
    else
        sendmail(mailto,subject,message);
    end
else
    sendmail(mailto,subject);
end
