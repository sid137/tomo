%ラゲール多項式
%2006/1/21 引数をベクトルで取れるようにした。

function LnX = Ln(N,alpha,x)

siz=size(x);
if N==0,
     LnX=1*ones(siz(1),siz(2));
     return;
else

  if N==1,
        LnX=1+alpha-x;
        return;
  end
  
end  

if N >1,

  if siz(1)==1 

    arrayLn=zeros(N+1,siz(2));
    arrayLn(1,:)=1;
    arrayLn(2,:)=1+alpha-x;
      for n=3:N+1,
        arrayLn(n,:)=((2*(n-1)-1+alpha-x).*arrayLn(n-1,:)...
                     -(n-2+alpha).*arrayLn(n-2,:))./(n-1);
      end
	LnX=arrayLn(N+1,:);
  else

   if siz(2)==1 

    arrayLn=zeros(siz(1),N+1);
    arrayLn(:,1)=1;
    arrayLn(:,2)=1+alpha-x;
      for n=3:N+1,
        arrayLn(:,n)=((2*(n-1)-1+alpha-x).*arrayLn(:,n-1)...
                     -(n-2+alpha).*arrayLn(:,n-2))./(n-1);
      end
	LnX=arrayLn(:,N+1);
   end
        if siz(1)~=1 & siz(2)~=1
	 disp(['引数xはベクトルのみ計算できます。'])
	return;
	end
  end

end


