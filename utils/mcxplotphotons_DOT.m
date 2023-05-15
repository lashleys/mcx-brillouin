function varargout=mcxplotphotons_DOT(traj,detpt,varargin)
%
%    mcxplotphotons(traj)
%       or
%    mcxplotphotons(traj, 'color','r','marker','o')
%    [sorted, linehandle]=mcxplotphotons(traj)
%
%    Plot photon trajectories from MCXLAB's output
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        traj: the 5th output of mcxlab, storing the photon trajectory info
%           traj.id: the photon index being recorded
%           traj.pos: the 3D position of the photon; for each photon, the
%                  positions are stored in serial order
%           traj.data: the combined output, in the form of 
%                [id,pos,weight,reserved]'
%
%    output:
%        sorted: a structure to store the sorted trajectory info
%        sorted.id: the sorted vector of photon id, staring from 0
%        sorted.pos: the sorted position vector of each photon, only
%                 recording the scattering sites.
%
%    this file is part of Monte Carlo eXtreme (MCX)
%    License: GPLv3, see http://mcx.sf.net for details
%

if(~isstruct(traj) && size(traj,2)==6)
    traj=struct('id',typecast(single(traj(:,1)),'uint32'),'pos',traj(:,2:4),'weight',traj(:,5));
end

[newid, idx]=sort(traj.id);
newpos=traj.pos(idx,:);
hg=plot3(newpos(:,1),newpos(:,2),newpos(:,3),'-',varargin{:});

det_pos = detpt.p;
[r, ~] = size(det_pos);
r = 4;
for i = 1:r
    x_ind = find(newpos(:,1) == det_pos(i,1));
    y_ind = find(newpos(:,2) == det_pos(i,2));
    z_ind = find(newpos(:,3) == det_pos(i,3));
    % Find common index
    res = intersect(x_ind, y_ind);
    ind = intersect(res,z_ind);
    % plot all points for the id that the index corresponds to
    id = newid(ind);
    det_newpos_x = newpos(newid == id, 1);
    det_newpos_y = newpos(newid == id, 2);
    det_newpos_z = newpos(newid == id, 3);
    plot3(det_newpos_x, det_newpos_y, det_newpos_z)
    hold on
    plot3(det_newpos_x(1), det_newpos_y(1), det_newpos_z(1), 'ro')
    plot3(det_newpos_x(end), det_newpos_y(end), det_newpos_z(end), 'go')
end

output={struct('id',newid, 'pos',newpos), hg};
[varargout{1:nargout}]=output{1:nargout};
