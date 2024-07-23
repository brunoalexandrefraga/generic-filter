n = 14;
fc = 2e3;

[zb,pb,kb] = butter(n,2*pi*fc,"s");
[bb,ab] = zp2tf(zb,pb,kb);
[hb,wb] = freqs(bb,ab,4096);

[ze,pe,ke] = ellip(n,3,30,2*pi*fc,"s");
[be,ae] = zp2tf(ze,pe,ke);
[he,we] = freqs(be,ae,4096);

plot([wb we]/(2e3*pi), ...
    mag2db(abs([hb he])))
axis([0 5 -60 5])
grid
xlabel("Frequency (GHz)")
ylabel("Attenuation (dB)")
legend(["butter" "ellip"])