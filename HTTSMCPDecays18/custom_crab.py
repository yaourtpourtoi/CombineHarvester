from CRABClient.UserUtilities import getUsernameFromCRIC
def custom_crab(config):
    print('>> Customising the crab config')
    config.General.workArea='Impacts-19052020'
    config.Site.storageSite = 'T2_UK_London_IC'

    config.JobType.maxMemoryMB = 3000

    config.Data.outLFNDirBase='/store/user/{}/{}/'.format(
        getUsernameFromCRIC(), config.General.workArea
    )
    config.JobType.allowUndistributedCMSSW = True
