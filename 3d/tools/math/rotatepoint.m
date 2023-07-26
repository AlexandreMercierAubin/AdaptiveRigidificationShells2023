function newp = rotatepoint(q,p)
    pQuat = quaternion([zeros(size(p,1),1),p]);
    qQuat = quaternion(q);
    newqTotal = qQuat * pQuat * conj(qQuat);
    newp = compact(newqTotal); 
end

