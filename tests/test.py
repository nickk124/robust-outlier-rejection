# for Travis
import rcr

data = [0, 0.1, -0.1, 0.2, -0.2, 1e50]

r = rcr.RCR(rcr.LS_MODE_68)
r.performBulkRejection(data) # perform outlier rejection

# View results
rejected_data = r.result.rejectedY

assert len(rejected_data) >= 1