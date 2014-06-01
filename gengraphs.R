library(ggplot2)

df = read.table('empty.txt', header=TRUE)

hyp = function(a,b,c,d) log((c-a)*(d-b)/(b-a)/(d-c))

pdf.a4r = function(file, ...)
    pdf(file, paper='a4r', width=16.54, height=11.69, ...)


pdf.a4r('hyperbolic.pdf')
print(qplot(c, hyp(vconj, p, q, u), data=df, geom='line',
            main="Hyperbolic length of return interval inside critical values"))
dev.off()

pdf.a4r('relative_length_L_R.pdf')
dg = data.frame(c=rep(df$c, 2),
                ratio=c((df$q-df$c)/df$u, (df$c-df$p)/(1-df$vconj)),
                side=rep(c('R/u', 'L/v'), each=nrow(df)))
print(qplot(c, ratio, data=dg, geom='line', col=side,
            main="Relative length of R in [0,u] and L in [vconj,1]"))

dg$ratio = c( (df$c-df$p)/(df$c-df$xi1), (df$q-df$c)/(df$eta1-df$c) )
dg$side = rep(c('L/[xi1,c]', 'R/[eta1,c]'), each=nrow(df))
print(qplot(c, ratio, data=dg, geom='line', col=side,
            main="Relative length of L in [xi1,c] and R in [c,eta1]"))
dev.off()

pdf.a4r('distortion.pdf')
dg = data.frame(c=rep(df$c, 2),
                distortion=c(df$dist_left, df$dist_right),
                branch=rep(c('left','right'), each=nrow(df)))
print(qplot(c, distortion, data=dg, geom='line', col=branch,
            main="Distortion of left and right branch of first-return map"))
dev.off()

pdf.a4r('crenorm.pdf')
print(qplot(c, (c-p)/(q-p), data=df, geom='line',
            main="Critical value of the renormalization"))
dev.off()
