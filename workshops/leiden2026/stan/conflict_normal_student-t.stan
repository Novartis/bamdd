data {
  real df1;
  real location1;
  real scale1;
  real df2;
  real location2;
  real scale2;
}
parameters {
  real mu1;
  real mu2;
  real mu;
}
model {
  target += student_t_lpdf(mu1 | df1, location1, scale1);
  target += student_t_lpdf(mu2 | df2, location2, scale2);
  target += student_t_lpdf(mu | df1, location1, scale1);
  target += student_t_lpdf(mu | df2, location2, scale2);
}
