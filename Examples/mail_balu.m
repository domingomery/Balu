test_mail = 1; % <- 1 para probar, y 0 para enviar el mail definitivo.

heads = {...
    'Estimado Nombre 1:',...
    'Estimado Nombre 2:',...
    'Estimado Nombre 3:',...
    'Estimado Nombre 4:',...
    'Estimado Nombre 5:',...
    'Estimado Nombre 6:',...
    };


if test_mail==0
    
    mails = {...
        'nombre_1@gmail.com',...
        'nombre_2@yahoo.com',...
        'nombre_3@gmail.com',...
        'nombre_4@hotmail.com',...
        'nombre_5@gmail.com',...
        'nombre_6@hotmail.com',...
};
    
else % direcciones de email de prueba, coloca aqui distintos mails tuyos
    mails = {...
        'my_mail_1@gmail.com',...
        'my_mail_2@gmail.com',...
        'my_mail_3@gmail.com',...
};
    
end

% NOTA: En cada fila un alumno, en cada columna una nota:


NOTA = [
5.0	7.0	4.0	4.0	4.0	5.0	
5.0	4.0	4.0	3.0	4.0	5.0	
5.0	4.0	3.0	5.0	4.0	5.0	
2.0	3.0	4.0	4.0	4.0	5.0	
5.5	4.0	4.0	4.5	5.0	4.0	
5.0	4.0	4.0	4.0	4.0	5.0	
];


dsg = [
    '1. Nota 1                        '
    '2. Nota 2                        '
    '3. Nota 3                        '
    '4. Nota 4                        '
    '5. Nota 5                        '
    '6. Nota Final                    '
    ];

p = [0.1 0.1 0.15 0.25 0.25 0.1 0.05 1 1]'*100; % porcentajes de cada columna


m = size(NOTA,2);

n = length(mails);

msg0 = {'Adjunto te estoy enviando el archivo corregido de ....', ' ', 'El resumen de notas es el siguiente:'};

my_gmail    = 'my_gmail@gmail.com';   % <- poner aqui tu gmail
my_password = 'hola123';              % <- poner aqui tu password
subject     = 'MY SUBJECT';           % <- poner aqui asunto
signature   = {'Saludos',' ','My Name',' ','----------------------------------','Prof. XXXXXXX','Departamento de XXXXXX','Universidad Catolica de Chile','http://xxxxxxxxxxxxxx','----------------------------------'};


for i=1:n
    s = num2fixstr(i,2);
    archivo_pdf = ['Archivo_' s '.pdf'];
    msj = zeros(m,48);
    for j=1:m
            msj(j,:) = sprintf('%s: %4.1f (%5.1f%%)',dsg(j,:),NOTA(i,j),p(j));
    end
    msg = char(msj);
    
    head_msg = heads{i};
    
    Bio_sendmail(my_gmail,...
        pass_word,...
        mails{i},...
        subject,...
        [head_msg ' ' msg0 ' ' msg(1,:) ' ' msg(2,:) ' ' msg(3,:) ' ' msg(4,:) ' ' msg(5,:) ' ' msg(6,:) ' ' signature],...
        archivo_pdf);
    % si no hay archivo adjunto, usar el comando Bio_sendmail sin el ultimo argumento
    
end

