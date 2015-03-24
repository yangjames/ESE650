clear all
close all

pos = [];
ts = [];
t0 = [];
for i = 1:171
    name = ['joints1/' 'joint' num2str(i) '.mat'];
    dat = load(name);
    pos = [pos; dat.pos];
    if i == 1
        ts = [ts; dat.ts];
        t0 = dat.t0;
    else
        ts = [ts; dat.ts + ts(end)];
    end
end

save('joints1.mat','pos','ts','t0')