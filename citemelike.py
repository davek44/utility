#!/usr/bin/env python
from optparse import OptionParser
import urllib2
import pdb, re, math, random

################################################################################
# citemelike.py
#
# Choose a random paper from my citeulike account to read.
################################################################################

user = 'dakelley'
star_factor = 10

url_re = re.compile('href="(/user/dakelley/article/\d+)"')
star_re = re.compile('src="/static/img/star(\d).png"')

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    parser.add_option('--tag', dest='tag', help='Choose a paper with the given tag')
    (options,args) = parser.parse_args()

    if options.tag:
        citeulike_url = 'http://www.citeulike.org/user/%s/tag/%s/order/to_read' % (user,options.tag)
    else:
        citeulike_url = 'http://www.citeulike.org/user/%s/order/to_read' % user

    # get papers
    papers = get_papers(citeulike_url)

    if len(papers) == 0:
        parser.error('No papers with the tag %s' % options.tag)

    # re-score stars
    papers = [(math.pow(stars/5.0,star_factor),paper) for (stars,paper) in papers]

    # choose random paper
    max_rand = sum([score for (score,paper) in papers])
    rand_score = random.uniform(0,max_rand)
    rand_tmp = 0.0
    for (score,paper) in papers:
        rand_tmp += score
        if rand_tmp > rand_score:
            print 'http://www.citeulike.org%s' % paper
            break
        
################################################################################
# get_papers
#
# Get all papers from the following base url
################################################################################
def get_papers(citeulike_url):
    papers = []

    page_num = 1
    no_read = True
    unread_found = True
    while no_read and unread_found:
        unread_found = False

        #f = urllib2.urlopen('%s/page/%d' % (citeulike_url,page_num))
        req = urllib2.Request('%s/page/%d' % (citeulike_url,page_num), headers={'User-Agent':"Magic Broswer"})
        f = urllib2.urlopen(req)
        cul_read = f.read()
        cul_text = ''.join(cul_read)
        cul_lines = cul_text.split('\n')

        for line in cul_lines:
            if line.find('class="title"') != -1:
                url_match = url_re.search(line)
                paper_url = url_match.group(1)
                unread_found = True

            elif line.find('/static/img/star') != -1:
                star_match = star_re.search(line)
                stars = int(star_match.group(1))
                papers.append((stars,paper_url))
            
            elif line.find('radio') == -1 and line.find('already read') != -1:
                papers = papers[:-1]
                unread = False
                break

        page_num += 1

    return papers

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
