function yout = warburgShort(v,w)
yout = zeros(length(w),2);  % allocated yout


objfcn = @(v) v(1)+v(2).*tanh(v(3)*(j.*w).^v(4))./(v(3)*(j.*w).^v(4));


yout(:,1) = real(objfcn(v));
yout(:,2) = imag(objfcn(v));

end

